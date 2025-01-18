import re 
import os
import shutil
from rdkit import Chem
from pyteomics import mgf


def extract_collision_energy(params):
    """
    Extracts the collision energy from the compound name if it exists.

    @param params: Dictionary of spectrum parameters
    @return: Collision energy as a string if found, otherwise None.
    """

    #if the collision energy is already specified
    if params.get("collision_energy") != None:
        if params["collision_energy"] != "":
            return params["collision_energy"]

    #we check if the collision energy is in the compound name
    compound_name = params.get('compound_name').lower()
    match = re.search(r'collisionenergy[:\s]*([\d]+)', compound_name, re.IGNORECASE)
    if match: return match.group(1)
    else: return None

def precurseur_correction(precursor):
    """
    Format the precursor name without charge and without special characters.
    
    @param precursor: the precursor name
    @return: a cleaned precursor name string without brackets, charges, or special characters
    """
    #remove the brackets [ ]
    precursor = re.sub(r'[\[\]]', '', precursor)
    #remove the charge 
    #precursor = re.sub(r'(\S+)\d*([+-])$', r'\1', precursor)
    #precursor = re.sub(r'(\S+)\d*([+-])$', r'\1', precursor)
    precursor = re.sub(r'(\S+)\d+[+-]$', r'\1', precursor)
    precursor = re.sub(r'(\S+)([+-])$', r'\1', precursor)
    #precursor = re.sub(r'(\S+)(\d+[+-])$', r'\1', precursor)
    return precursor


def canonisation_smiles(smiles):
    """
    Canonizes a SMILES string using RDKit.

    @param smiles: SMILES string.
    @return: Canonized SMILES string if valid, otherwise None.
    @raises: Exception if RDKit fails to process the SMILES string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            return None
    except Exception as e:
        print(f"Error canonizing SMILES {smiles}: {e}")
        return None


def write_in_file(spectre, file, params):
    """
    Writes the corrected spectrum into the file in MGF format, corresponding to its collision energy and precursor.

    @param spectre: Spectrum data, including m/z and intensity arrays.
    @param file: File to write the spectrum to.
    @param params: Dictionary of spectrum parameters
    """

    file.write("BEGIN IONS\n")

    #we rewrite in the file everything that is not in a field to be removed
    for key, value in params.items():
        if key not in fields_to_remove:
            if isinstance(value, (tuple, list)):
                if value[1] != None:
                    value = ' '.join(map(str, value))
                else:
                    #prevents writing the 'None' value that appears in the pepmass field
                    value = value[0]

            file.write(f"{key.upper()}={value}\n")
    
    #we rewrite all the peaks/intensities of the spectrum
    mz_value = spectre['m/z array']
    intensity_value = spectre['intensity array']
    for mz, intensity in zip(mz_value, intensity_value):
        file.write(f"{mz} {intensity}\n")
    
    file.write("END IONS\n\n")


def clear_directory(directory):
    """
    Deletes the contents of the directory to prevent overwriting the molecules in the file during a new execution of the code.
    
    @param directory: Path to the directory to be cleared
    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)


def parse_mgf(filename, newdir, dico_smiles=None):
    """
    Parses an MGF file, processes spectra, and writes them into categorized files.
    
    @param filename: Path to the MGF file to parse
    @param newdir: Directory where the categorized MGF files will be saved
    @param dico_smiles: Dictionary to store SMILES strings categorized by file name
    
    @return: Updated dictionary containing SMILES strings grouped by file name
    """

    if not os.path.exists(newdir):
        os.mkdir(newdir)

    #clear_directory(newdir)
    if dico_smiles == None:
        dico_smiles = dict()

    with mgf.MGF(filename) as spectres:
        for spectre in spectres:
            #for each spectrum, we check the precursor, the SMILES, the collision energy, and its validity, then we write it in the corresponding file
            params = spectre['params']
            collision_energy = extract_collision_energy(params)
            precursortype = params.get(precursorname)
            precursortype = precurseur_correction(precursortype)
            smiles = params.get('smiles')
            smiles_canonise = canonisation_smiles(smiles)
            
            if "not validated" in params['compound_name']:
                newfilename = "not_validated.mgf"

            else:
                if (collision_energy is not None) and (precursortype != "") and smiles != "" and smiles_canonise != None:
                    params['collision_energy'] = collision_energy
                    params['compound_name'] = re.sub(r'collisionenergy[:\s]*([\d]+)', '', params['compound_name'], flags=re.IGNORECASE).strip()
                    newfilename = "energy_" + collision_energy + "_precursor_" + str(precursortype) + ".mgf"
                    dico_smiles[newfilename] = dico_smiles.get(newfilename, [])+[smiles_canonise]
            
                elif (collision_energy is None and smiles != ""):
                    newfilename = "no_energy.mgf"
            
                elif (smiles == "" and collision_energy is not None):
                    newfilename = "no_smiles.mgf"
            
                else:
                    #(no energy and no smiles) or no precursor 
                    newfilename = "inexploitable.mgf"
            
            newpath = newdir + os.sep + newfilename

            with open(newpath, 'a') as newfile:
                write_in_file(spectre, newfile, params)


    return dico_smiles


def write_smiles_to_file(dico_smiles, smiles_file):
    """
    Writes the SMILES dictionary to a text file.

    @param dico_smiles: Dictionary where keys are file names and values are lists of SMILES strings
    @param smiles_file: Path to the text file where SMILES will be written
    """

    with open(smiles_file, 'w') as sf:
        for key, value in dico_smiles.items():
            sf.write(key + '\n')
            for smiles in value:
                sf.write(smiles + '\n')




if __name__ == "__main__":

    filename = "Cluster.mgf"
    precursorname = "precursortype"
    newdir = "cluster_molecules_test"
    fields_to_remove = ["charge", "ontology", "ionmode", "instrumenttype", "instrument", "manufacturer", "ms_mass_analyzer", "ms_ionisation", "ms_dissociation_method", "scans"]
    # the precursor fields in cluster.mgf and all_gnps do not have the same name
    if filename == "Cluster.mgf" or filename == "mol.mgf": precursorname = "precursortype"
    else: precursorname = "adduct"
    dico_smiles = parse_mgf(filename, newdir)
    filename = "ALL_GNPS_cleaned.mgf"
    precursorname = "adduct"
    parse_mgf(filename, newdir, dico_smiles)
    path_smilesfile = newdir + os.sep + "smiles.txt"

    write_smiles_to_file(dico_smiles, path_smilesfile)
