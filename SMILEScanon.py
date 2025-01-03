from rdkit import Chem
import os

dirname = "molecules_parsees"
filename = ["no_energy.mgf", "no_smiles.mgf"]

#attention y a pas tout le temps l'inchikey (avec pyteomics si ce champs n'existe pas
#donne cette valeur "inchikey inconnu i")
#rajouter no smiles et no validated dans filename
#utiliser withopen

#la je crée juste un dictionnaire avec les SMILES canonisés associé à leur INCHKEY, mais
#mais peut-être qu'il faudra écrire le nouveau SMILE dans les fichiers .mgf
#et aussi le code est brouillon parce que je galère trop à apprendre python :/

def parse_into_dictionnary(dirname, filename):
    inch_smiles = {}
    
    for f in os.listdir(dirname):
        if f not in filenames :
            path = dirname + "/" + f
            with open(path, 'r') as file :
                content = file.readlines()
            inch = None
            smiles = None
            #ce serait plus joli avec un regexp mais je suis vraiment un caca en python
            for l in content:
                #on ne verifie pas que les champs ne sont pas vides, pt il faut mettre des and
                if l[0:8] == "INCHIKEY":
                    inch = l[9:].replace("\n", "")
                elif l[0:6] == "SMILES":
                    smiles = l[7:].replace("\n","")
                    inch_smiles[inch] = smiles
        
    return inch_smiles

def canonisation(inch_smiles):
    compt_diff_smiles = 0
    for inch, smiles in inch_smiles.items():
        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            inch_smiles[inch] = "SMILES invalide"
        else:
            tmp_smiles = Chem.MolToSmiles(mol, canonical=True)
            if tmp_smiles != smiles:
                #print(f"SMILES canonique : {tmp_smiles}")
                #print(f"SMILES du fichier : {smiles}")
                compt_diff_smiles += 1
            inch_smiles[inch] = tmp_smiles
    return compt_diff_smiles

def verif_with_inchkey(inch_smiles):
    compt_smiles_invalides = 0
    compt_diff_inch = 0
    for inch, smiles in inch_smiles.items():
        if smiles == "SMILES invalide":
            compt_smiles_invalides += 1
        else:
            mol = Chem.MolFromSmiles(smiles)
            inch_from_smile = Chem.inchi.MolToInchiKey(mol)
            if inch_from_smile != inch:
                compt_diff_inch += 1
    return compt_smiles_invalides, compt_diff_inch

inch_smiles = parse_into_dictionnary(dirname, filename)
compt_diff_smiles = canonisation(inch_smiles)
compt_smiles_invalides, compt_diff_inch = verif_with_inchkey(inch_smiles)
print(f"nombre de smiles différents : {compt_diff_smiles}")
print(f"nombre de smiles invalides : {compt_smiles_invalides}")
print(f"nombre d'Inchikey différentes : {compt_diff_inch}")
