import re 
import os
import shutil
from pyteomics import mgf

filename = "Cluster.mgf"
newdir = "molecules_parsees"

def extract_collision_energy(params):
    """
    Extracts the collision energy from the compound name if it exists.
    """

    compound_name = params.get('compound_name').lower()
    match = re.search(r'collisionenergy[:\s]*([\d]+)', compound_name, re.IGNORECASE)
    if match: return match.group(1)
    else: return None


def write_in_file(spectre, file, params):
    """
    Writes the corrected spectrum into the file corresponding to its collision energy and precursor.
    """

    file.write("BEGIN IONS\n")

    for key, value in params.items():
        if isinstance(value, (tuple, list)):  # Cas de PEPMASS
            value = ' '.join(map(str, value))
        file.write(f"{key.upper()}={value}\n")

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

    clear_directory(newdir)

    with mgf.MGF(filename) as spectres:
        for spectre in spectres:
            params = spectre['params']
            collision_energy = extract_collision_energy(params)
            precursortype = params.get('precursortype')


            if (collision_energy is not None):
                params['collision_energy'] = collision_energy
                params['compound_name'] = re.sub(r'collisionenergy[:\s]*([\d]+)', '', params['compound_name'], flags=re.IGNORECASE).strip()
                newfilename = "energy_" + collision_energy + "_precursor_" + str(precursortype) + ".mgf"
            else:
                newfilename = "no_energy.mgf"
            
            newpath = newdir + os.sep + newfilename

            with open(newpath, 'a') as newfile:
                write_in_file(spectre, newfile, params)


parse_mgf(filename, newdir)
