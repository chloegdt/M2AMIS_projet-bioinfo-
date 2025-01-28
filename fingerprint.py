import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from parse_data_final import create_dict_from_smiles


FILENAME = "molecules_parsees/energy_102040_precursor_M+H.mgf"
FPSIZE = 4096


def getMolecules(filename):
    """ Récupère les smiles du fichier et les transformes en molécules """
    smiles = list()
    with mgf.MGF(filename) as reader:
        for molecule in reader:
            smiles.append(molecule['params']['smiles'])
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def getMorganFingerprintsFromFile(filename):
    """ Calcule les fingerprints de Morgan des molecules du fichier """
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    return molecules, [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]

def getMorganFingerprints(smiles):
    """ Calcule les fingerprints de Morgan à partir d'une liste de SMILES """
    molecules = [Chem.MolFromSmiles(smile) for smile in smiles]
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    return [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]

def getRDkitFingerprintsFromFile(filename):
    """ Calcule les fingerprints RDkit des molecules du fichier """
    rdkitGen = rdFingerprintGenerator.GetRDKitFPGenerator(minPath=1, maxPath=7, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    return molecules, [np.array(rdkitGen.GetFingerprint(molecule)) for molecule in molecules]


def fingerprintsSimilarity(fingerprints):
    """
        Calcule la similarité de Tanimoto entre chaque fingerprints
        Tanimoto(A, B) = (A et B) / (A ou B)
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
    """ Prend une liste de liste (matrice) en entrée et la sauvegarde au format csv dans le fichier similarity/nom.csv """
    folder = "similarity/"
    os.makedirs(folder, exist_ok=True)
    file = folder + name + '.csv'

    array = np.array(matrix)
    np.savetxt(file, array, fmt='%.6f', delimiter=',')
    return file


def matrixToTxt(matrix, name):
    folder = "data/"
    os.makedirs(folder, exist_ok=True)
    file = folder + name + '.txt'
    with open(file, 'w') as f:
        for i, list_sim in enumerate(matrix):
            for j, sim in enumerate(list_sim):
                if sim > 0:
                    f.write(f"{i+1} {j+1} {sim}\n")


def createEveryMatrix(filename):
    """ Créé toutes les matrices de similarité """
    smiles_dict = create_dict_from_smiles(filename)
    for filename, smiles in smiles_dict.items():
        fg = getMorganFingerprints(smiles)
        fg_sim = fingerprintsSimilarity(fg)
        matrixToCSV(fg_sim, f"fgsim_{filename}")



if __name__ == '__main__':
    molecules, fg = getMorganFingerprintsFromFile("cluster_molecules_test/energy_25.0_precursor_M+H.mgf")
    # molecules, fg = getRDkitFingerprints(FILENAME)
    fg_sim = fingerprintsSimilarity(fg)
    matrixToTxt(fg_sim, "fingerprints_energy_25.0_precursor_M+H.mgf")
    print(fg_sim)

    # print(sorted(fg_sim[-1],reverse=False))
    # createEveryMatrix("cluster_molecules_test/smiles.txt")
