import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
import hdbscan
import matplotlib.pyplot as plt
import umap

class SpectraClustering:
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


    def perform_clustering(self):
        dense_matrix = self.matrix.toarray()
        reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
        umap_embedding = reducer.fit_transform(self.distance_matrix)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=5, metric='manhattan', alpha=1.0)
        self.labels = clusterer.fit_predict(umap_embedding)

        print("Matrice de distance (10 premières lignes) :")
        print(self.distance_matrix[:10, :10])
        print("Labels de clusters :")
        print(np.unique(self.labels, return_counts=True))
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        print(f"Nombre de clusters détectés : {n_clusters}")
        return self.labels


    def visualize_clusters(self):
        reducer = umap.UMAP(n_components=5, metric='precomputed', random_state=42)
        umap_results = reducer.fit_transform(self.distance_matrix)

        plt.figure(figsize=(10, 8))
        plt.scatter(umap_results[:, 0], umap_results[:, 1], c=self.labels, cmap='viridis', s=50, alpha=0.7)
        plt.title('Clustering des Spectres avec HDBSCAN et UMAP')
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


file = 'resultats_spectres_cosinus/energy_50.0_precursor_M+H.txt'

clustering = SpectraClustering(file)
clustering.load_data()
labels = clustering.perform_clustering()
clustering.save_clusters_to_file('energy_50.0_precursor_M+H.txt')

print("Labels de clusters :")
print(labels)

#clustering.visualize_clusters()

i = 0
for l in labels:
    if l == -1:
        i += 1

print(i)