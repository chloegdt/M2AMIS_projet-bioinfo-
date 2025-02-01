import numpy as np
from scipy.sparse import coo_matrix
import markov_clustering as mc
import matplotlib.pyplot as plt
import random
import os
import shutil

#YOU NEED TO INSTALL markov_clustering BEFORE RUNNING THE PROGRAM
#EXAMPLE : pip install markov_clustering[drawing]

def create_matrix(file):
    rows = []
    cols = []
    data = []

    with open(file, 'r') as f:
        for line in f:
            try:
                spectre1, spectre2, similarity = line.split()
                spectre1, spectre2, similarity = int(spectre1) - 1, int(spectre2) - 1, float(similarity)

                rows.append(spectre1)
                cols.append(spectre2)
                data.append(similarity)
            except Exception as e:
                continue

        matrix = coo_matrix((data, (rows, cols)))
    return matrix

def clear_directory(directory):
    """
    Deletes the contents of the directory to prevent overwriting the molecules in the file during a new execution of the code.
    
    @param directory: Path to the directory to be cleared
    """
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

def draw_cluster(outputdir, file, matrix, clusters):
    outputdir_img = outputdir + "_graph"
    if not os.path.exists(outputdir_img):
        os.mkdir(outputdir_img)
        
    output_file_img = os.path.join(outputdir_img, f"{os.path.splitext(file)[0]}_graph.png")
    mc.draw_graph(matrix, clusters, node_size=50, with_labels=False, edge_color="silver")
    plt.savefig(output_file_img, format='png')  # Sauvegarder en PNG
    plt.close()

def write_cluster(clusters, file):
    with open(file, 'w') as f:
        for cluster in clusters:
            for value in cluster:
                f.write(str(value + 1) + " ")
            f.write("\n")

def clustering(inputdir, outputdir, inflation, draw_graph):
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
        matrix = create_matrix(os.path.join(inputdir, file))
        result = mc.run_mcl(matrix, inflation = inflation)
        clusters = mc.get_clusters(result)
        output_file = os.path.join(outputdir, f"{os.path.splitext(file)[0]}_cluster.txt")
        write_cluster(clusters, output_file)

        if draw_graph:
            draw_cluster(outputdir, file, matrix, clusters) 

if __name__ == "__main__":
    inputdir = "resultats_spectres"
    inflation = 2.0
    outputdir = "resultats_clusters"
    clustering(inputdir, outputdir, inflation, True)