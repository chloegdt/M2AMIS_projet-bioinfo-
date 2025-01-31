import subprocess
import os
import shutil

#YOU NEED TO INSTALL MCL BEFORE RUNNING THE PROGRAM
#EXAMPLE : (Linux) sudo apt-get install mcl

def clear_directory(directory):
    """
    Deletes the contents of the directory to prevent overwriting the molecules in the file during a new execution of the code.
    
    @param directory: Path to the directory to be cleared
    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)


def clustering(inputdir, outputdir, inflation):
    """
    Calls the command mcl to create clusters (using markov clustering) from inputdir.

    @param inputdir: Path to the directory containing similarity matrix
    @param outputdir: Path to the directory containing results (for each file, each line is a cluster)
    @param influation: when increasing, the max distance between two elements of a cluster is smaller
    """
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    
    clear_directory(outputdir)

    for file in os.listdir(inputdir):
        commande = ["mcl", os.path.join(inputdir, file), "--abc", "-I", inflation, "-odir", outputdir]
        resulta  = subprocess.run(commande)

if __name__ == "__main__":
    # inputdir = "resultats_smiles_tanimoto"
    inputdir = "resultats_spectres_cosinus"
    # outputdir = "clusters_smiles_tanimoto"
    outputdir = "clusters_spectres_cosinus"
    # inflation = "14.0" # pour les smiles
    # inflation = "1.1" # pour les spectres
    inflation = "1.1"
    clustering(inputdir, outputdir, inflation)