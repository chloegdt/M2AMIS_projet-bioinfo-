import os
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from parse_data_final import create_dict_from_smiles
from fingerprint import fingerprintSimilarity, matrixToTxt

FILENAME = "cluster_molecules_test/smiles.txt"


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


def getFunctionalGroupCounts(molecules):
    """
    Calcule le nombre d'occurrences de chaque groupe fonctionnel pour chaque molécule.
    
    @param molecules: liste d'objets rdkit molécule.
    @return: liste de dictionnaires de groupes fonctionnels en clés et le nombre d'apparition en valeur.
    """
    functional_groups = []
    catalog = FilterCatalog.GetFunctionalGroupHierarchy()

    for mol in molecules:
        groups = {}
        matches = catalog.GetMatches(mol)

        for match in matches:
            group_name = match.GetDescription()
            groups[group_name] = groups.get(group_name, 0) + 1
        functional_groups.append(groups)
    return functional_groups


def getGroupMatrix(functional_groups):
    """
    Transforme chaque dictionnaire en liste pour pouvoir calculer la similarité.
    
    @param functional_groups: liste de dictionnaires de groupes fonctionnels en clés et le nombre d'apparition en valeur.
    @return: matrice contenant les groupes fonctionnels de chaque molécule.
    """
    # On récupère la liste de tous les groupes fonctionnels contenus dans les molécules de l'entrée.
    groups_name = set()
    for groups in functional_groups:
        groups_name.update(groups.keys())
    groups_name = sorted(groups_name)

    # Matrice (nombre de molécules x nombre de groupes fonctionnels utilisés)
    matrix = np.zeros((len(functional_groups), len(groups_name)))
    for i, groups in enumerate(functional_groups):
        for j, group_name in enumerate(groups_name):
            matrix[i, j] = groups.get(group_name, 0)

    return matrix, groups_name


def manhattanDistance(functional_groups):
    """
    Calcule la matrice de distances Manhattan entre les molécules à partir de leurs groupes fonctionnels.
    
    @param functional_groups: matrice contenant les groupes fonctionnels de chaque molécule.
    @return: matrice des distances Manhattan.
    """
    n = len(functional_groups)
    distance_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            distance = np.sum(np.abs(functional_groups[i] - functional_groups[j]))
            distance_matrix[i, j] = distance
    return distance_matrix






def TEMPgetFunctionalGroups(molecules):
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

    molecules_matches = np.zeros((len(molecules), len(functional_groups)), dtype=int)

    for i, mol in enumerate(molecules):
        atoms_used = set()
        for j, group in enumerate(functional_groups):
            matches = mol.GetSubstructMatches(group)
            for match in matches:
                if not atoms_used.isdisjoint(match):
                    continue
                molecules_matches[i][j] = 1
                atoms_used.update(match)

    return molecules_matches





if __name__ == "__main__":
    f1 = "cluster_molecules_test/energy_0.0_precursor_M+H.mgf"
    f2 = "cluster_molecules_test/energy_10.0_precursor_M+H.mgf"
    file1 = "cluster_molecules_test/energy_50.0_precursor_M+H.mgf"
    file2 = "cluster_molecules_test/energy_10.0_precursor_M+H.mgf"
    file3 = "cluster_molecules_test/energy_30.0_precursor_M+H.mgf"
    # mol = getMolecules(f2)
    # func = getFunctionalGroupCounts(mol)
    # matrix, groups = getGroupMatrix(func)
    # dist = manhattanDistance(matrix)
    # print(dist)
    # print(func)
    # print(groups)

    """
    filters = FilterCatalog.GetFlattenedFunctionalGroupHierarchy()
    for k, pat in filters.items():
        print(f"{k} -> {Chem.MolToSmarts(pat)}")
    """



    molecules = getMolecules(file1)
    groupes = TEMPgetFunctionalGroups(molecules)
    sim = fingerprintSimilarity(groupes)
    matrixToTxt(sim, "cluster_molecules/resultats_groupes_fonc/", file1.split('/')[1])
    print("Fichier 1")

    molecules = getMolecules(file2)
    groupes = TEMPgetFunctionalGroups(molecules)
    sim = fingerprintSimilarity(groupes)
    matrixToTxt(sim, "cluster_molecules/resultats_groupes_fonc/", file2.split('/')[1])
    print("Fichier 2")

    molecules = getMolecules(file3)
    groupes = TEMPgetFunctionalGroups(molecules)
    sim = fingerprintSimilarity(groupes)
    matrixToTxt(sim, "cluster_molecules/resultats_groupes_fonc/", file3.split('/')[1])
    print("Fichier 3")
