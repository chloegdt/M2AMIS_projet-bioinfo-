import os
from matchms.importing import load_from_mgf
import matchms.filtering as ms_filters
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matchms import calculate_scores
from matchms.Spectrum import Spectrum
from matchms.similarity.BaseSimilarity import BaseSimilarity


class ManhattanSimilarity(BaseSimilarity):
    """
    Fonction de calcul de la similarité Manhattan entre deux spectres.
    """
    # Définir le type de données attendu pour les scores
    score_datatype = np.dtype([('score', 'f4'), ('matched_peaks', 'i4')])

    def pair(self, reference: Spectrum, query: Spectrum):
        # similarité avec la même molécule
        if reference is query:
            return np.array((1.0, len(np.union1d(reference.peaks.mz, query.peaks.mz))), dtype=self.score_datatype)
        
        # Aligner les spectres par m/z
        ref_mz, ref_intensities = reference.peaks.mz, reference.peaks.intensities
        query_mz, query_intensities = query.peaks.mz, query.peaks.intensities

        # Créer un vecteur commun basé sur les m/z présents dans les deux spectres
        common_mz = np.union1d(ref_mz, query_mz)

        # Interpolation des intensités sur les m/z alignés
        ref_interp = np.interp(common_mz, ref_mz, ref_intensities, left=0, right=0)
        query_interp = np.interp(common_mz, query_mz, query_intensities, left=0, right=0)

        # Calcul de la distance Manhattan
        manhattan_distance = np.sum(np.abs(ref_interp - query_interp))
        
        # Calcul de la similarité inversée (valeur principale entre 0 et 1)
        similarity = 1 / (1 + manhattan_distance)
        
        # Calcul du nombre de pics communs
        matched_peaks = len(common_mz)
        
        # Retourner un tableau structuré compatible avec matchms
        return np.array((similarity, matched_peaks), dtype=self.score_datatype)


def metadata_processing(spectrum):
    """
    Filtrage et standardization des spectre.
    """
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.repair_inchi_inchikey_smiles(spectrum)
    spectrum = ms_filters.harmonize_undefined_smiles(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchi(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchikey(spectrum)
    spectrum = ms_filters.add_precursor_mz(spectrum)
    return spectrum


def peak_processing(spectrum):
    """
    Filtrage et normalisation des pics.
    """
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.normalize_intensities(spectrum)
    spectrum = ms_filters.select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = ms_filters.select_by_mz(spectrum, mz_from=10, mz_to=1000)
    return spectrum


def save_to_csv(scores, output_csv_path):
    """
    Sauvegarde dans un fichier CSV.
    """
    data = []
    for (reference, query, score) in scores:
        data.append({"reference": reference.metadata.get("title", reference.metadata.get("formula", "unknown")), # spectre principale
                     "query": query.metadata.get("title", query.metadata.get("formula", "unknown")),  # spectre secondaire 
                     "score": score, # score de la similarité et nombre de pics qui match
                     })

    df_scores = pd.DataFrame(data)
    df_scores.to_csv(output_csv_path, index=False)

    print(f"Manhattan similarity scores saved to {output_csv_path}")


####### MAIN ########
if __name__ == "__main__":
    input_directory = "cluster_molecules_test"  # Répertoire contenant les fichiers .mgf
    output_directory = "similarite_results"      # Répertoire pour sauvegarder les résultats

    # Fichiers à ignorer
    files_to_ignore = {"inexploitable.mgf", "no_energy.mgf", "no_smiles.mgf", "not_validated.mgf", "smiles.txt"}

    # Créer le répertoire de sortie s'il n'existe pas
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Parcourir tous les fichiers .mgf dans le répertoire
    for mgf_file in os.listdir(input_directory):
        if mgf_file.endswith(".mgf") and mgf_file not in files_to_ignore:
            input_file_path = os.path.join(input_directory, mgf_file)
            output_csv_path = os.path.join(output_directory, f"{os.path.splitext(mgf_file)[0]}_manhattan_scores.csv")

            print(f"Processing file: {mgf_file}")

            # Charger les spectres depuis le fichier .mgf
            spectrums = list(load_from_mgf(input_file_path))

            # Appliquer les traitements aux spectres
            spectrums = [metadata_processing(s) for s in spectrums]
            spectrums = [peak_processing(s) for s in spectrums]

            # Calculer la similarité Manhattan
            similarity_measure = ManhattanSimilarity()
            scores = calculate_scores(spectrums, spectrums, similarity_measure, is_symmetric=True)

            # Sauvegarder les scores dans un fichier CSV
            save_to_csv(scores, output_csv_path)

        else:
            print(f"Skipping file: {mgf_file}")


# # test sur un seul fichier 
# if __name__ == "__main__":
#     input_directory = "cluster_molecules_test"  # Répertoire contenant les fichiers .mgf
#     output_directory = "test_manhattan_results"  # Répertoire pour sauvegarder les résultats

#     selected_file = "energy_0.0_precursor_M-H.mgf"
#     input_file_path = os.path.join(input_directory, selected_file)

#     if not os.path.exists(input_file_path):
#         print(f"Erreur : le fichier '{selected_file}' n'existe pas dans le répertoire '{input_directory}'.")
#     else:
#         if not os.path.exists(output_directory):
#             os.makedirs(output_directory)

#         output_csv_path = os.path.join(output_directory, f"{os.path.splitext(selected_file)[0]}_manhattan_scores.csv")
#         print(f"Processing file: {selected_file}")

#         spectrums = list(load_from_mgf(input_file_path))
#         spectrums = [metadata_processing(s) for s in spectrums]
#         spectrums = [peak_processing(s) for s in spectrums]

#         similarity_measure = ManhattanSimilarity()
#         scores = calculate_scores(spectrums, spectrums, similarity_measure, is_symmetric=True)

#         save_to_csv(scores, output_csv_path)
#         print(f"Résultats sauvegardés dans : {output_csv_path}")
