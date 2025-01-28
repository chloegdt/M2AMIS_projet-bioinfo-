import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from parse_data_final import create_dict_from_smiles

FILENAME = "cluster_molecules_test/smiles.txt"
FPSIZE = 4096


def getSmiles(filename):
    """
    Récupère les smiles d'un fichier mgf.
    
    @param filename: chemin du fichier mgf à lire.
    @return: liste de smiles.
    """
    smiles = list()
    with mgf.MGF(filename) as reader:
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


def fingerprintSimilarity(fingerprints):
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
            fp_or = np.sum(np.bitwise_or(fingerprints[i], fingerprints[j]))
            fp_and = np.sum(np.bitwise_and(fingerprints[i], fingerprints[j]))

            similarity = (fp_and / fp_or) if (fp_or) else 0
            similarity_matrix[i, j] = similarity
    return similarity_matrix


def matrixToTxt(matrix, filename):
    """
    Sauvegarde une matrice dans un fichier .txt sous forme de matrice creuse.

    @param matrix: matrice de similarité.
    @param filename: nom du fichier.
    """
    directory_path = "cluster_molecules/resultats_fingerprints/"
    output_path = os.path.join(directory_path, f"{os.path.splitext(os.path.basename(filename))[0]}.txt")
    os.makedirs(directory_path, exist_ok=True)

    with open(output_path, 'w') as f:
        for i, list_sim in enumerate(matrix):
            for j, sim in enumerate(list_sim):
                if sim > 0:
                    f.write(f"{i+1} {j+1} {sim}\n")


def createEveryMatrix(file_path):
    """
    Calcule et sauvegarde (en fichier txt) les matrices de similarité à partir d'un fichier .txt contenant les smiles.
    
    @param filename: chemin du fichier .txt à lire.
    """
    smiles_dict = create_dict_from_smiles(file_path)

    for filename, smiles in smiles_dict.items():
        fingerprints = getMorganFingerprints(smiles)
        matrix = fingerprintSimilarity(fingerprints)
        matrixToTxt(matrix, f"fingerprints_{filename.split('.')[0]}")



if __name__ == '__main__':
    # file_test = "cluster_molecules_test/energy_100.0_precursor_M+H.mgf"

    # fingerprints = getMorganFingerprints(getSmiles(file_test))
    # matrix = fingerprintSimilarity(fingerprints)
    # matrixToTxt(matrix, "fingerprints_" + file_test.split("/")[1].split(".")[0])

    createEveryMatrix(FILENAME)
