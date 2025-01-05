import re 
import os
import shutil
from pyteomics import mgf


def extract_collision_energy(params):
    """
    Extracts the collision energy from the compound name if it exists.
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
    """
    #remove the brackets [ ]
    precursor = re.sub(r'[\[\]]', '', precursor)
    #remove the charge 
    precursor = re.sub(r'(\S+)(\d+[+-])$', r'\1', precursor)
    return precursor

def write_in_file(spectre, file, params):
    """
    Writes the corrected spectrum into the file corresponding to its collision energy and precursor.
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
    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)


def parse_mgf(filename, newdir):
    """
    Main function for parsing the data and rewriting it into new files, correcting the compound name, and filling in the collision energy field.
    """

    if not os.path.exists(newdir):
        os.mkdir(newdir)

    #clear_directory(newdir)

    with mgf.MGF(filename) as spectres:
        for spectre in spectres:
            #for each spectrum, we check the precursor, the SMILES, the collision energy, and its validity, then we write it in the corresponding file
            params = spectre['params']
            collision_energy = extract_collision_energy(params)
            precursortype = params.get(precursorname)
            precursortype = precurseur_correction(precursortype)
            smiles = params.get('smiles')
            
            if "not validated" in params['compound_name']:
                newfilename = "not_validated.mgf"

            else:
                if (collision_energy is not None) and (precursortype != "") and smiles != "":
                    params['collision_energy'] = collision_energy
                    params['compound_name'] = re.sub(r'collisionenergy[:\s]*([\d]+)', '', params['compound_name'], flags=re.IGNORECASE).strip()
                    newfilename = "energy_" + collision_energy + "_precursor_" + str(precursortype) + ".mgf"
            
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



filename = "Cluster.mgf"
#filename = "ALL_GNPS_cleaned.mgf"
newdir = "cluster_molecules4"
fields_to_remove = ["charge", "ontology", "ion_mode", "instrumenttype", "instrument", "ms_manufacturer", "ms_mass_analyzer", "ms_ionisation", "ms_dissociation_method", "scans"]
# the precursor fields in cluster.mgf and all_gnps do not have the same name
if filename == "Cluster.mgf" or filename == "mol.mgf": precursorname = "precursortype"
else: precursorname = "adduct"


parse_mgf(filename, newdir)
