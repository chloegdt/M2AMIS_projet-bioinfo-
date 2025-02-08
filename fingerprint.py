import os
import logging
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from parse_data_final import create_dict_from_smiles

FILENAME = "cluster_molecules/smiles.txt"
SAVE_DIRECTORY = "cluster_molecules/resultats_fingerprints/"
FPSIZE = 4096
FILES = [
    "energy_50.0_precursor_M+H.mgf",
    "energy_10.0_precursor_M+H.mgf",
    "energy_30.0_precursor_M+H.mgf",
    "energy_25.0_precursor_M+Na.mgf",
    "energy_37.0_precursor_M+Na.mgf",
]


def getSmiles(file_path):
    """
    Récupère les smiles d'un fichier mgf.
    
    @param file_path: chemin du fichier mgf à lire.
    @return: liste de smiles.
    """
    smiles = list()
    with mgf.MGF(file_path) as reader:
        for molecule in reader:
            smiles.append(molecule['params']['smiles'])
    return smiles


def getMorganFingerprints(smiles):
    """
    Calcule les fingerprints de Morgan à partir d'une liste de SMILES.
    
    @param smiles: liste de smiles.
    @return: liste de fingerprints (numpy array binaire).
    """
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles]
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    return [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]


def tanimotoSimilarity(fingerprints):
    """
    Calcule la similarité de Tanimoto entre chaque fingerprints.
    Tanimoto(A, B) = (A et B) / (A ou B)
    
    @param fingerprints: liste de fingerprints (numpy array binaire).
    @return: matrice de similarité.
    """
    n = len(fingerprints)
    similarity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if i == j:
                similarity_matrix[i, j] = 1
                continue

            fp_or = np.sum(np.bitwise_or(fingerprints[i], fingerprints[j]))
            fp_and = np.sum(np.bitwise_and(fingerprints[i], fingerprints[j]))

            similarity = (fp_and / fp_or) if (fp_or) else 0
            similarity_matrix[i, j] = similarity
    return similarity_matrix


def matrixToTxt(matrix, directory_path, filename):
    """
    Sauvegarde une matrice de similarité dans un fichier .txt sous forme de matrice creuse.

    @param matrix: matrice de similarité.
    @param directory_path: chemin du dossier de sauvegarde.
    @param filename: nom du fichier de sauvegarde.
    """
    output_path = os.path.join(directory_path, f"{os.path.splitext(os.path.basename(filename))[0]}.txt")
    os.makedirs(directory_path, exist_ok=True)

    with open(output_path, 'w') as f:
        for i, list_sim in enumerate(matrix):
            for j, sim in enumerate(list_sim):
                if sim > 0:
                    f.write(f"{i+1} {j+1} {sim}\n")


def createEveryMatrix(file_path, directory_path):
    """
    Calcule et sauvegarde (en fichier txt) les matrices de similarité à partir d'un fichier .txt contenant les smiles.
    
    @param file_path: chemin du fichier .txt à lire.
    @param directory_path: chemin du dossier de sauvegarde.
    """
    logging.info(f"Début du traitement de tous les fichiers.")
    smiles_dict = create_dict_from_smiles(file_path)

    for filename, smiles in smiles_dict.items():
        fingerprints = getMorganFingerprints(smiles)
        matrix = tanimotoSimilarity(fingerprints)
        matrixToTxt(matrix, directory_path, filename)

    logging.info(f"Création des matrices de similarité terminé (résultats dans le dossier: {directory_path})")


def main(files):
    """
    Fonction principale pour calculer et sauvegarder (en fichier txt) les matrices de similarité des fichiers spécifiés en entrée.

    @param files: Liste contenant les nom des fichiers à traiter.
    """
    if not os.path.exists(FILENAME):
        logging.error(f"Fichier {FILENAME} introuvable. Executez le parsing puis réessayez.")
        return

    smiles_dict = create_dict_from_smiles(FILENAME)
    smiles = [smiles_dict.get(file, []) for file in files]
    logging.info(f"Début du calcul des similarités fingerprints des SMILES.")

    for i, smile in enumerate(smiles):
        fingerprints = getMorganFingerprints(smile)
        matrix = tanimotoSimilarity(fingerprints)
        matrixToTxt(matrix, SAVE_DIRECTORY, files[i])
        logging.info(f"Fichier {files[i]} traité.")

    logging.info(f"Calcul des similarités fingerprints des SMILES terminé. \nRésultats dans le dossier: {SAVE_DIRECTORY}")



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    main(FILES)
    # createEveryMatrix(FILENAME, SAVE_DIRECTORY)
