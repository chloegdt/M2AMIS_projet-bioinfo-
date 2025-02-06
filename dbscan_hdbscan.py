import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors
import hdbscan
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import umap

class DensityClustering:
    def __init__(self, file):
        self.file = file
        self.matrix = None
        self.labels = None

    def load_data(self):
        rows = []
        cols = []
        data = []

        with open(self.file, 'r') as f:
            for line in f:
                spectre1, spectre2, similarity = line.split()
                spectre1, spectre2, similarity = int(spectre1) - 1, int(spectre2) - 1, float(similarity)

                rows.append(spectre1)
                cols.append(spectre2)
                data.append(similarity)

        self.matrix = coo_matrix((data, (rows, cols)))
        self.matrix = self.matrix + self.matrix.T - coo_matrix((self.matrix.data, (self.matrix.col, self.matrix.row)), shape=self.matrix.shape)
        self.distance_matrix = 1 - self.matrix.toarray()

    def find_epsilon(self, min_samples):

        neigh = NearestNeighbors(n_neighbors=min_samples)
        nbrs = neigh.fit(self.distance_matrix)
        distances, indices = nbrs.kneighbors(self.distance_matrix)
        distances = np.sort(distances, axis=0)
        distances = distances[:,1]
        plt.plot(distances)

        plt.title('Valeur optimale d\'epsilon pour DBSCAN')
        plt.show()

    def perform_clustering(self, clustering):
        dense_matrix = self.matrix.toarray()
        reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
        umap_embedding = reducer.fit_transform(self.distance_matrix)
        if clustering == "HDBSCAN":
            clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=4, metric='euclidean', alpha=1.0)
        elif clustering == "DBSCAN":
            clusterer = DBSCAN(eps=0.2, min_samples=4, metric='euclidean')
        else:
            print("ERREUR : ni DBSCAN ni HDBSCAN")
            return self.labels
        
        self.labels = clusterer.fit_predict(umap_embedding)

        print("Matrice de distance (10 premières lignes) :")
        print(self.distance_matrix[:10, :10])
        print("Labels de clusters :")
        print(np.unique(self.labels, return_counts=True))
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        print(f"Nombre de clusters détectés : {n_clusters}")
        return self.labels


    def visualize_clusters(self, clustering):
        reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
        umap_results = reducer.fit_transform(self.distance_matrix)

        plt.figure(figsize=(10, 8))
        plt.scatter(umap_results[:, 0], umap_results[:, 1], c=self.labels, cmap='viridis', s=50, alpha=0.7)
        plt.title(f'Clustering des Spectres avec {clustering} et UMAP')
        plt.xlabel('UMAP Dimension 1')
        plt.ylabel('UMAP Dimension 2')
        plt.colorbar(label='Cluster ID')
        plt.show()

    def save_clusters_to_file(self, output_file):
        clusters_dict = {}
        for idx, label in enumerate(self.labels):
            if label != -1:
                if label not in clusters_dict:
                    clusters_dict[label] = []
                clusters_dict[label].append(idx + 1)

        with open(output_file, 'w') as f:
            for cluster_id, ids in clusters_dict.items():
                f.write("\t".join(map(str, ids)) + "\n")


if __name__ == '__main__':
    file = 'resultats_spectres_cosinus/energy_50.0_precursor_M+H.txt'

    clustering = DensityClustering(file)
    clustering.load_data()
    #clustering.find_epsilon(4)
    labels = clustering.perform_clustering("DBSCAN")
    clustering.save_clusters_to_file('clusters_sepctres/energy_50.0_precursor_M+H.txt')

    print("Labels de clusters :")
    print(labels)
    print("Bruit :")
    print(sum(1 for l in labels if l == -1))

    clustering.visualize_clusters("DBSCAN")