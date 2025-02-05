import os
import logging
import numpy as np

from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from parse_data_final import create_dict_from_smiles
from fingerprint import tanimotoSimilarity, matrixToTxt, createEveryMatrix, FILENAME, create_dict_from_smiles


def getFunctionalGroups(molecules):
    """
    Calcule une matrice contenant pour chaque molécule la présence d'une liste de groupes fonctionnels (1 présent, 0 non).

    @param molecules: liste d'objets molécules rdkit.
    @return: matrice contenant les groupes fonctionnels de chaque molécule.
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


def createEveryMatrix(file_path, directory_path):
    """
    Calcule et sauvegarde (en fichier txt) les matrices de similarité à partir d'un fichier .txt contenant les smiles.

    @param file_path: chemin du fichier .txt à lire.
    @param directory_path: chemin du dossier de sauvegarde.
    """
    logging.info(f"Début du traitement de tous les fichiers.")
    smiles_dict = create_dict_from_smiles(file_path)

    for filename, smiles in smiles_dict.items():
        molecules = [Chem.MolFromSmiles(smile) for smile in smiles]
        groupes = getFunctionalGroups(molecules)
        matrix = tanimotoSimilarity(groupes)
        matrixToTxt(matrix, directory_path, filename)
    logging.info(f"Création des matrices de similarité terminé (résultats dans le dossier: {directory_path})")



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    files = [
        "energy_50.0_precursor_M+H.mgf",
        "energy_10.0_precursor_M+H.mgf",
        "energy_30.0_precursor_M+H.mgf",
        "energy_25.0_precursor_M+Na.mgf",
    ]

    smiles_dict = create_dict_from_smiles(FILENAME)
    smiles = [smiles_dict.get(file) for file in files]

    for i, smile in enumerate(smiles):
        molecules = [Chem.MolFromSmiles(s) for s in smile]
        groupes = getFunctionalGroups(molecules)
        sim = tanimotoSimilarity(groupes)
        matrixToTxt(sim, "cluster_molecules/resultats_groupes-fonc/", files[i])
        logging.info(f"Fichier {i+1} traité.")
    logging.info(f"Traitement terminé, résultats dans le dossier: cluster_molecules/resultats_groupes-fonc/")

    # createEveryMatrix(FILENAME, "cluster_molecules/resultats_groupes-fonc/")
