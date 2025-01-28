import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from parse_data_final import create_dict_from_smiles


FILENAME = "cluster_molecules_test/smiles.txt"
FPSIZE = 4096


def getMolecules(filename):
    """
    Récupère les smiles d'un fichier mgf et les transformes en objets rdkit molécule.
    
    @param filename: chemin du fichier mgf à lire.
    @return: liste d'objets molécule.
    """
    smiles = list()
    with mgf.MGF(filename) as reader:
        for molecule in reader:
            smiles.append(molecule['params']['smiles'])
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def getMorganFingerprintsFromFile(filename):
    """
    Calcule les fingerprints de Morgan des molecules à partir d'un fichier mgf.
    
    @param filename: chemin du fichier mgf à lire.
    @return molecules: liste d'objets molécule.
    @return fingerprints: liste de fingerprints (numpy array binaire).
    """
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    return molecules, [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]

def getMorganFingerprints(smiles):
    """
    Calcule les fingerprints de Morgan à partir d'une liste de SMILES.
    
    @param smiles: liste de smiles.
    @return : liste de fingerprints (numpy array binaire).
    """
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles]
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    return [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]

def getRDkitFingerprintsFromFile(filename):
    """
    Calcule les fingerprints RDkit des molecules à partir d'un fichier mgf.
    
    @param filename: chemin du fichier mgf à lire.
    @return molecules: liste d'objets molécule.
    @return fingerprints: liste de fingerprints (numpy array binaire).
    """
    rdkitGen = rdFingerprintGenerator.GetRDKitFPGenerator(minPath=1, maxPath=7, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    return molecules, [np.array(rdkitGen.GetFingerprint(molecule)) for molecule in molecules]


def fingerprintsSimilarity(fingerprints):
    """
    Calcule la similarité de Tanimoto entre chaque fingerprints.
    Tanimoto(A, B) = (A et B) / (A ou B)
    
    @param fingerprints: liste de fingerprints (numpy array binaire).
    @return: matrice de similarité.
    """
    n_fp = len(fingerprints)
    similarity_matrix = np.zeros((n_fp, n_fp))

    for fp in range(n_fp):
        for target_fp in range(fp, n_fp):
            fp_or = np.sum(np.bitwise_or(fingerprints[fp], fingerprints[target_fp]))
            fp_and = np.sum(np.bitwise_and(fingerprints[fp], fingerprints[target_fp]))

            similarity = (fp_and / fp_or) if (fp_or) else 0
            # similarity_matrix[fp, target_fp] = similarity
            similarity_matrix[target_fp, fp] = similarity
    return similarity_matrix


def matrixToCSV(matrix, name):
    """
    Sauvegarde une matrice dans un fichier .csv sous forme de matrice pleine.

    @param matrix: matrice de similarité.
    @param name: nom du fichier.
    """
    folder = "similarity/"
    os.makedirs(folder, exist_ok=True)
    file = folder + name + '.csv'

    array = np.array(matrix)
    np.savetxt(file, array, fmt='%.6f', delimiter=',')


def matrixToTxt(matrix, name):
    """
    Sauvegarde une matrice dans un fichier .txt sous forme de matrice creuse.

    @param matrix: matrice de similarité.
    @param name: nom du fichier.
    """
    # TODO: changer nom folder et utiliser os.join path
    folder = "data/"
    os.makedirs(folder, exist_ok=True)
    file = folder + name + '.txt'
    with open(file, 'w') as f:
        for i, list_sim in enumerate(matrix):
            for j, sim in enumerate(list_sim):
                if sim > 0:
                    f.write(f"{i+1} {j+1} {sim}\n")


def createEveryMatrix(filename):
    """
    Calcule et sauvegarde (en fichier CSV) les matrices de similarité à partir d'un fichier .txt.
    
    @param filename: chemin du fichier .txt à lire.
    """
    smiles_dict = create_dict_from_smiles(filename)
    for filename, smiles in smiles_dict.items():
        fg = getMorganFingerprints(smiles)
        fg_sim = fingerprintsSimilarity(fg)
        matrixToCSV(fg_sim, f"fgsim_{filename}")



if __name__ == '__main__':
    file_test = "cluster_molecules_test/energy_100.0_precursor_M+H.mgf"

    molecules, fg = getMorganFingerprintsFromFile(file_test)
    fg_sim = fingerprintsSimilarity(fg)
    matrixToTxt(fg_sim, "fingerprints_" + file_test.split("/")[1])

    # createEveryMatrix(FILENAME)
