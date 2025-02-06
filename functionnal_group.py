import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from parse_data_final import create_dict_from_smiles
from fingerprint import tanimotoSimilarity, matrixToTxt, createEveryMatrix, FILENAME


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


def getFunctionalGroups(molecules):
    """
    Calcule la matrice de distances Manhattan entre les molécules à partir de leurs groupes fonctionnels.
    
    @param molecules: liste d'objet molécules rdkit.
    @return: matrice des distances Manhattan.
    """
    functional_groups = [
        "[CX3](=O)[OX2H1]", # Acide Carboxylique
        "[CX3](=O)O[C](=O)[#6]", # Anhydride d'acide
        "[CX3](=O)[Cl,Br,I]", # Halogénure d'acyle
        "[NX3][CX3](=[OX1])[#6]", # Amide
        "[CX3](=O)[OX2H0][#6]", # Ester
        "[CX3](=O)[#6]", # Cétone
        "[CX3H1](=O)[#6]", # Aldéhyde
        "[OX2H]", # Alcool
        "[OD2]([#6])[#6]", # Ether
        "c1ccccc1", # Arène / Cycle aromatique
        "[CX4][F,Cl,Br,I]", # Halogénoalcane
        "[NX3;H2]", #Amine
        "[#6]-[#6]", # Alcane
        "C#C", # Alcyne
        "C=C", #Alcène
    ]
    functional_groups = [Chem.MolFromSmarts(smarts) for smarts in functional_groups]

    molecules_groups = np.zeros((len(molecules), len(functional_groups)), dtype=int)

    for i, mol in enumerate(molecules):
        atoms_used = set()
        for j, group in enumerate(functional_groups):
            matches = mol.GetSubstructMatches(group)
            for match in matches:
                if not atoms_used.isdisjoint(match):
                    continue
                molecules_groups[i][j] = 1
                atoms_used.update(match)

    return molecules_groups


def manhattanDistance(molecules_groups):
    """
    Calcule la distance Manhattan entre les molécules à partir de leurs groupes fonctionnels.
    
    @param molecules_groups: matrice contenant les groupes fonctionnels de chaque molécule.
    @return: matrice des distances Manhattan.
    """
    n = len(molecules_groups)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            distance = np.sum(np.abs(molecules_groups[i] - molecules_groups[j]))
            distance_matrix[i, j] = distance
    return distance_matrix



if __name__ == "__main__":
    files = [
        "cluster_molecules_test/energy_50.0_precursor_M+H.mgf",
        "cluster_molecules_test/energy_10.0_precursor_M+H.mgf",
        "cluster_molecules_test/energy_30.0_precursor_M+H.mgf",
    ]

    for i, file in enumerate(files):
        molecules = getMolecules(file)
        groupes = getFunctionalGroups(molecules)
        sim = tanimotoSimilarity(groupes)
        matrixToTxt(sim, "cluster_molecules/resultats_groupes-fonc/", file.split('/')[1])
        print(f"Fichier {i+1} traité.")


    # createEveryMatrix(FILENAME, "cluster_molecules/resultats_groupes-func/")
