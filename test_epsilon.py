import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

def load_data(file):
        rows = []
        cols = []
        data = []

        with open(file, 'r') as f:
            for line in f:
                spectre1, spectre2, similarity = line.split()
                spectre1, spectre2, similarity = int(spectre1) - 1, int(spectre2) - 1, float(similarity)

                rows.append(spectre1)
                cols.append(spectre2)
                data.append(similarity)

        matrix = coo_matrix((data, (rows, cols)))
        matrix = matrix + matrix.T - coo_matrix((matrix.data, (matrix.col, matrix.row)), shape=matrix.shape)
        distance_matrix = 1 - matrix.toarray()
        return distance_matrix

def find_epsilon(distance_matrix, min_samples):

    neigh = NearestNeighbors(n_neighbors=min_samples)
    nbrs = neigh.fit(distance_matrix)
    distances, indices = nbrs.kneighbors(distance_matrix)
    distances = np.sort(distances, axis=0)
    distances = distances[:,1]
    plt.plot(distances)

    plt.title('Valeur optimale d\'epsilon pour DBSCAN')
    plt.show()

matrix = load_data("resultats_smiles_tanimoto/fg_energy_50.0_precursor_M+H.txt")
find_epsilon(matrix, 8192)