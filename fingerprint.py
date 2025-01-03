import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


FILENAME = "molecules_parsees/energy_102040_precursor_M+H.mgf"
FPSIZE = 4096


def getMolecules(filename):
    """ Récupère les smiles du fichier et les transformes en molécules """
    smiles = set()
    with mgf.MGF(filename) as reader:
        for molecule in reader:
            smiles.add(molecule['params']['smiles'])
    return [Chem.MolFromSmiles(smile) for smile in smiles]

def getMorganFingerprints(filename):
    """ Calcule les fingerprints de Morgan des molecules du fichier """
    morganGen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    # molecules.append(Chem.MolFromSmiles("c1c(cccc1)CCO"))
    # molecules.append(Chem.MolFromSmiles("c1c(cccc1)CCCO"))
    return molecules, [np.array(morganGen.GetFingerprint(molecule)) for molecule in molecules]

def getRDkitFingerprints(filename):
    """ Calcule les fingerprints RDkit des molecules du fichier """
    rdkitGen = rdFingerprintGenerator.GetRDKitFPGenerator(minPath=1, maxPath=7, fpSize=FPSIZE)
    molecules = getMolecules(filename)
    # molecules.append(Chem.MolFromSmiles("c1c(cccc1)CCO"))
    # molecules.append(Chem.MolFromSmiles("c1c(cccc1)CCCO"))
    return molecules, [np.array(rdkitGen.GetFingerprint(molecule)) for molecule in molecules]


def fingerprintsSimilarity(fingerprints):
    """
    Calcule la similarité de Tanimoto entre chaque fingerprints
    Tanimoto(A, B) = (A et B) / (A ou B)
    """
    n_fp = len(fingerprints)
    similarity_matrix = np.zeros((n_fp, n_fp))

    for fp in range(n_fp):
        for target_fp in range(fp + 1, n_fp):
            fp_or = np.sum(np.bitwise_or(fingerprints[fp], fingerprints[target_fp]))
            fp_and = np.sum(np.bitwise_and(fingerprints[fp], fingerprints[target_fp]))

            similarity = (fp_and / fp_or) if (fp_or) else 0
            similarity_matrix[fp, target_fp] = similarity
            similarity_matrix[target_fp, fp] = similarity
    return similarity_matrix



if __name__ == '__main__':
    # molecules, fg = getMorganFingerprints(FILENAME)
    molecules, fg = getRDkitFingerprints(FILENAME)
    fg_sim = fingerprintsSimilarity(fg)

    print(fg_sim)
    # print(sorted(fg_sim[-1],reverse=False))
