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
        commande = ["mcl", os.path.join(inputdir, file), "--abc", "-I", influation, "-odir", outputdir]
        resulta  = subprocess.run(commande)

if __name__ == "__main__":
    inputdir = "resultats_spectres"
    inflation = "2.0"
    outputdir = "resultats_clusters"
    clustering(inputdir, outputdir, influation)