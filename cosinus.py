import os
import multiprocessing
import logging
import numpy as np
import pandas as pd

from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from pathlib import Path
from matchms.importing import load_from_mgf
from matchms.similarity import CosineGreedy
from matchms import calculate_scores
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_intensity
from matchms import Spectrum

INPUT_DIRECTORY = "cluster_molecules/"
SAVE_DIRECTORY = "cluster_molecules/resultats_cosinus_spectres/"


def get_mz_values(spectres):
    """
    Récupère la plus petite valeur de m/z et la plus grande parmi tous les spectres du fichier

    @param spectres: Liste d'objets Spectrum de matchms
    @return min_mz, max_mz: la plus petite valeur et la plus grande valeur de m/z
    """
    # initialise le min et le max à 0
    min_mz = 0
    max_mz = 0

    for spectre in spectres:
        if spectre.peaks:
            mz_values = spectre.peaks.mz
            min_mz = min(min_mz, mz_values.min())
            max_mz = max(max_mz, mz_values.max())

    # prévient l'utilisateur si aucune valeur n'a été trouvée
    if min_mz == 0 or max_mz == 0:
        logging.warning("Aucune valeur de m/z valide n'a été trouvée")

    return min_mz, max_mz


def extract_compoundname(spectre):
    """
    Extrait le nom du composé chimique depuis les métadonnées d'un spectre.

    @param spectre: objet Spectrum de matchms
    @return: nom du composé ou "Unknown" si absent
    """
    compound_name = spectre.metadata.get("compound_name", "Unknown")
    if compound_name == "Unknown":
        logging.warning("Nom de composé manquant pour un spectre.")
    return compound_name
    

def assign_ids_to_spectra(spectres):
    """
    Assigne un identifiant unique à chaque spectre dans la liste.
    
    @param spectres: Liste d'objets Spectrum de matchms
    @return: Liste de spectres avec des IDs uniques dans leur métadonnées
    """
    # si la liste est vide, un message d'erreur s'affiche
    if not spectres:
        logging.warning("La liste des spectres est vide. Aucun ID assigné.")
        return

    # ajoute un id unique dans les métadonnées du spectre
    for idx, spectre in enumerate(spectres):
        spectre.set("id", idx + 1)
        # logging.info(f"ID assigné : {spectre.metadata.get('id')}")
    return spectres


def process(spectre, step=1, threshold=0.1):
    """
    Filtre et normalise un spectre en utilisant les intensités.

    @param spectre: le spectre à traiter
    @return spectre : le spectre filtré et normalisé
    """
    
    # si il y a des pics, on peut agréger
    if spectre.peaks:
        mz_values = spectre.peaks.mz
        intensity_values = spectre.peaks.intensities
        
        # on définit des valeurs de m/z selon un pas
        rounded_mz = (step * (mz_values / step).round()).astype(float)

        # Agrégation des intensités en fonction des valeurs de m/z arrondies
        unique_mz = np.unique(rounded_mz)
        aggregated_intensities = [np.sum(intensity_values[rounded_mz == mz]) for mz in unique_mz]
            
        # Mise à jour des valeurs m/z et des intensités agrégées
        
        new_peaks = np.vstack([unique_mz, aggregated_intensities]).T
        spectre = Spectrum(
            mz=new_peaks[:, 0],
            intensities=new_peaks[:, 1],
            metadata=spectre.metadata
        )

        # normalisation des intensités 
        spectre = normalize_intensities(spectre)
        # garde uniquement les pics dont l'intensité est supérieure au seuil choisi
        spectre = select_by_intensity(spectre, intensity_from=threshold)
    
    # si il ne reste plus de pics conservés, un message s'affiche
    if not spectre.peaks.mz.size:
        logging.warning("Le spectre est vide après le filtrage d'intensité.")
        return None
    
    return spectre


def save_score_to_file(scores_data, output_file):
    """
    Enregistre les score dans un fichier csv

    @param scores_data: liste des scores à sauvegarder
    @param output_file: chemin du fichier de sortie 
    """

    directory = os.path.dirname(output_file)
    if not os.path.exists(directory):
        os.makedirs(directory,exist_ok=True)
        print(f"Répertoire {directory} créé.")
 
    # Extraire indices et scores à partir des dictionnaires
    rows = [item["Ref_id"] for item in scores_data]
    cols = [item["Query_id"] for item in scores_data]
    scores = [item["Similarité"] for item in scores_data]
    scores = [score if score<int(1) else 1 for score in scores]

    # Créer la matrice creuse
    sparse_matrix = coo_matrix((scores, (rows, cols)))

    # Sauvegarder en format txt
    with open(output_file, 'w') as f:
        for i, j, v in zip(sparse_matrix.row, sparse_matrix.col, sparse_matrix.data):
            f.write(f"{i+1} {j+1} {v}\n")

    logging.info(f"Score sauvegardés dans le fichier : {output_file}")


def compute_cosine_scores(spectres, tolerance):
    """
    Calcule les scores de la similarité cosinus entre les spectres.

    @param spectre: object Spectrum de matchms
    @param tolerance: tolérance utilisée dans le calcul
 
    @return score: l'object score de matchms
    @return n: nombre de spectres dans le fichier
    """
    n = len(spectres)
    cosine_greedy = CosineGreedy(tolerance=tolerance)
    score = calculate_scores(references=spectres, queries=spectres, similarity_function=cosine_greedy, is_symmetric=False)
    logging.info(f"Scores calculés pour {len(spectres)} spectres.")
    return n, score


def process_scores(n, score):
    """
    Extrait les informations des scores et structure les résultats.

    @param n: nombre de spectres dans le fichier
    @param score: score calculé par la fonction compute_cosine_scores
    @return scores_data: liste structurée des scores 
    """
    scores_data = [] 
    for query, ref, cosine_score in score:
        query_id = query.metadata.get("id", "Unknown") - 1 
        ref_id = ref.metadata.get("id", "Unknown") - 1

        if query_id <= ref_id:
            # logging.info(f"Spectres comparés : {query_id} et {ref_id}")
            scores_data.append({
                "Ref_id": ref_id,
                "Query_id": query_id,
                "Similarité":cosine_score[0]
                })
        
    return scores_data


def calculate(spectres, output_file="output.txt", tolerance=0.01):
    """
    Calcule les scores de similarité Cosinus entre les spectres et les sauvegardes dans un fichier.
    
    @param spectres: liste d'objets Spectrum de matchms 
    @param output_file: nom du fichier de sortie pour les scores
    @param tolerance: tolérance utilisée dans le calcul
    @return: liste des données de scores
    """

    # calcule des scores
    n, score = compute_cosine_scores(spectres, tolerance)
    
    # liste des scores
    scores_data = process_scores(n, score)

    # sauvegarde les scores dans un fichier 
    save_score_to_file(scores_data, output_file) 
    
    return scores_data


def process_file(file_path, directory_path, tolerance=0.01):
    """
    Calcule et sauvegarde la similarité cosinus entre les spectres d'un fichier.

    @param file_path: Fichier MGF des spectres à comparer.
    @param directory_path: Répertoire utilisé pour la sauvegarde.
    @param tolerance : Tolérance utilisée dans le calcul.
    """
    logging.info(f"Traitement du fichier : {file_path}")
    
    fichier = list(load_from_mgf(file_path))
    if not fichier:
        logging.error(f"Aucun spectre chargé depuis le fichier : {file_path}")
        return
        
    spectres = [process(spectre) for spectre in fichier if process(spectre) is not None]
    spectres = assign_ids_to_spectra(spectres)
    output_file = os.path.join(directory_path, "resultats_cosinus_spectres", f"{os.path.splitext(os.path.basename(file_path))[0]}.txt")
    
    calculate(spectres, output_file=output_file, tolerance=tolerance)



def main(directory_path, tolerance=0.01):
    """
    Fonction principale.

    @param file_path: fichier MGF des spectres à comparer
    @param output_file: chemin du fichier de sortie
    @param tolerance : tolérance utilisée dans le calcul
    """

    if not Path(directory_path).exists() or not Path(directory_path).is_dir():
        logging.error(f"Répertoire introuvable : {directory_path}")
        return

    target_directory = os.path.join(directory_path, "resultats_spectres")
    if not os.path.exists(target_directory):
        os.makedirs(target_directory, exist_ok=True)
        logging.info(f"Répertoire {target_directory} créé.")
    existing_files = [os.path.splitext(file)[0] for file in os.listdir(target_directory)]
    files = [f for f in os.listdir(directory_path)if f.startswith("energy") and f.endswith(".mgf") and os.path.splitext(f)[0] not in existing_files]
    
    if not files:
        logging.error("Aucun fichier n'a été trouvé'")
        return


    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        tasks = [(os.path.join(directory_path, file), directory_path, tolerance) for file in files]
        pool.starmap(process_file, tasks)
        


def main_selected_files(selected_files, directory_path=INPUT_DIRECTORY, tolerance=0.01):
    """
    Fonction principale pour traiter uniquement les fichiers spécifiés.
    
    @param directory_path: Répertoire contenant les fichiers à traiter.
    @param selected_files: Liste contenant les nom des fichiers à traiter.
    @param tolerance: Tolérance utilisée dans le calcul.
    """
    if not Path(directory_path).exists() or not Path(directory_path).is_dir():
        logging.error(f"Répertoire introuvable : {directory_path}")
        return

    files = [file for file in selected_files if os.path.exists(os.path.join(directory_path, file))]

    if not files:
        logging.error("Aucun fichier n'a été trouvé.")
        return

    logging.info(f"Début du calcul des similarités cosinus des spectres.")
    for file in files:
        process_file(os.path.join(directory_path, file), directory_path, tolerance)
    logging.info(f"Calcul des similarités cosinus des spectres terminé.")



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    directory_path = "cluster_molecules"
    main(directory_path)
