import numpy as np
import os
from parse_data_final import create_dict_from_smiles
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

def cluster_into_list(filename):
    """ Transform clusters from a file in a list of clusters 
    For example, file with 1\t5\t4\n2\t6\n3 is change into [[1,5,4][2,6][3]]

    @param filename: Path to file where clusters are written
    """
    cluster_list = []
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            parsing = line.split("\t")
            cluster_list.append([int(p) for p in parsing])
    return cluster_list

def id_to_smiles(filename, cluster_list, dico_smiles):
    """ Use the dictionnary from parse_data_final.py to create a list of cluster
    with SMILES instead of IDs

    @param filename: key of the dictionnary (format : energy_XX.X_precursor_X.mgf)
    @param cluster_list: list of cluster with IDs
    @param dico_smiles: dictionnary where the key is a filename from smiles.txt and the rest is 
    a list of SMILES
    """
    cluster_smiles_list =[]
    for clusters in cluster_list:
        smiles_list = []
        for id in clusters:
            smiles_list.append(dico_smiles[filename][id - 1])
        cluster_smiles_list.append(smiles_list)
    return cluster_smiles_list

def cluster_to_molecules(cluster_smiles_list):
    """ Use rdkit to transform a list of clusters with SMILES into a list of cluters
    with molecules

    @param cluster_smiles_list: list of clusters with SMILES
    """
    molecules_list = []
    for cluster in cluster_smiles_list:
        molecule_list = []
        for smiles in cluster:
            molecule_list.append(Chem.MolFromSmiles(smiles))
        molecules_list.append(molecule_list)
    return molecules_list

def getFunctionalGroups(molecules):
    """ Returns the number of occurences and the name of the functional group the most seen in a
    cluster

    @param molecules: a list of molecules in which we want to count the occurences of some functional groups
    """
    functional_groups_common_name = [
        "Acide Carboxylique", # Acide Carboxylique
        "Anhydride d'acide", # Anhydride d'acide
        "Halogénure d'acyle", # Halogénure d'acyle
        "Amide", # Amide
        "Ester", # Ester
        "Cétone", # Cétone
        "Aldéhyde", # Aldéhyde
        "Alcool", # Alcool
        "Ether", # Ether
        "Arène / Cycle aromatique", # Arène / Cycle aromatique
        "Halogénalcane", # Halogénoalcane
        "Amine", #Amine
        "Alcane", # Alcane
        "Alcyne", # Alcyne
        "Alcène", #Alcène
    ]
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

    ind = 0
    maxi = 0
    for col in range(len(functional_groups)):
        somme = 0
        for line in range(len(molecules)):
            somme += molecules_groups[line][col]
        if somme > maxi:
            maxi = somme
            ind = col
    return functional_groups_common_name[ind], int(maxi)

def groupCluster(molecules_list):
    """ For each cluster, return the functional group the most present, its number of occurences
    and the number of element in the cluster

    @param molecules_list: list of clusters with molecules
    """
    cluster_functionalGroup_list = []
    for molecules in molecules_list:
        best_functional_group, occurence = getFunctionalGroups(molecules)
        cluster_functionalGroup_list.append([len(molecules), best_functional_group, occurence])
    return cluster_functionalGroup_list
        

if __name__ == "__main__":
    file_smilesdico = os.path.join("cluster_molecules", "smiles.txt")
    dico_smiles = create_dict_from_smiles(file_smilesdico)
    file_cluster = os.path.join("cluster_molecules","resultats_clusters_hdbscan","groupefonct","energy_50.txt")
    cluster_list = cluster_into_list(file_cluster)
    cluster_smiles_list = id_to_smiles("energy_50.0_precursor_M+H.mgf", cluster_list, dico_smiles)
    molecules_list = cluster_to_molecules(cluster_smiles_list)
    print(groupCluster(molecules_list))