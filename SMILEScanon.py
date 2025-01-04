from pyteomics import mgf
from rdkit import Chem
import os

dirname = "cluster_molecules4"
filenames = ["no_energy.mgf", "no_smiles.mgf", "inexploitable.mgf", "not_validated.mgf"]

# Créer un dictionnaire avec les SMILES canonisés associés à leur INCHKEY
# Est-ce qu'on doit stocker quelque part le fichier associé ?
def parse_into_dictionnary(dirname, filenames):
    """
    Extrait un dictionnaire associant chaque InChIKey à son SMILES à partir d'un fichier MGF.
    Si l'InChIKey est manquant, il est nommé 'InChIKey inconnu X'.
    
    Args:
    - dirname (str): chemin du dossier contenant les molécules parsées
    - filenames (list): fichiers n'ayant pas de SMILES, ou pas d'énergie de collision
    
    Returns:
    - inch_smiles: un dictionnaire où les clés sont les InChIKeys et les valeurs les SMILES.
    """
    # inchi_counter = 1  # Compteur pour les InChIKeys manquants
    inch_smiles = {}
    
    for f in os.listdir(dirname):
        if f not in filenames:
            #path = os.path.join(dirname, test)
            path = os.path.join(dirname, f)
            with mgf.MGF(path) as spectres:
                for spectre in spectres:
                    params = spectre['params']
                    inch = params.get('inchikey')
                    smiles = params.get('smiles')

                    # Si l'InChIKey est manquant, attribuer un InChIKey inconnu X
                    # Alors les molécules identiques sont séparé par leur energie de collision + precurseur
                    if inch is None:
                        inch = f"InChIKey inconnu {params.get('compound_name')}"
                        # inch = f"InChIKey inconnu {inchi_counter}"
                        # inch += f"+ fichier : {f}" #Pour retrouver dans quel fichier est une molécule
                        # inchi_counter += 1

                    inch_smiles[inch] = smiles
        
    return inch_smiles

# Canoniser les SMILES
def canonisation(inch_smiles):
    """
    Canonise chaque SMILES de molécule et met à jour le dictionnaire avec les SMILES canonisés
    
    Args:
    - inch_smiles: un dictionnaire dont les SMILES ne sont pas canonisés
    
    Returns:
    - inch_smiles: un dictionnaire où les clés sont les InChIKeys et les valeurs les SMILES canonisés.
    - compt_diff_smiles: un compteur des SMILES qui n'étaient pas canonisés
    """
    compt_diff_smiles = 0
    for inch, smiles in inch_smiles.items():
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            inch_smiles[inch] = "SMILES invalide"
            print(smiles)
            print(inch)
        else:
            tmp_smiles = Chem.MolToSmiles(mol, canonical=True)
            if tmp_smiles != smiles:
                compt_diff_smiles += 1
            inch_smiles[inch] = tmp_smiles

    return inch_smiles, compt_diff_smiles

# Vérifier la cohérence des SMILES avec les InChIKeys
def verif_with_inchkey(inch_smiles):
    """
    Effectue des statistiques pour vérifier la cohérence des SMILES canonisés avec les InChIKeys.
    
    Args:
    - inch_smiles: un dictionnaire de SMILES canonisés
    
    Returns:
    - compt_smiles_invalides: compte les SMILES invalides
    - compt_diff_inchkeys: compte les InChIKeys qui ne correspondent pas au SMILES
    """
    compt_smiles_invalides = 0
    compt_diff_inchkeys = 0
    for inch, smiles in inch_smiles.items():
        if smiles == "SMILES invalide":
            compt_smiles_invalides += 1
        else:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                if not inch.startswith("InChIKey inconnu"):
                    inch_from_smile = Chem.MolToInchiKey(mol)
                    if inch_from_smile != inch:
                        compt_diff_inchkeys += 1
    return compt_smiles_invalides, compt_diff_inchkeys

inch_smiles = parse_into_dictionnary(dirname, filenames)
inch_smiles, compt_diff_smiles = canonisation(inch_smiles)
compt_smiles_invalides, compt_diff_inchkeys = verif_with_inchkey(inch_smiles)
print(f"Nombre de SMILES différents : {compt_diff_smiles}")
print(f"Nombre de SMILES invalides : {compt_smiles_invalides}")
print(f"Nombre d'InChIKeys différents : {compt_diff_inchkeys}")
with open('output.txt', 'w') as f:
    for inch, smiles in inch_smiles.items():
        f.write(f"{inch}: {smiles}\n")
