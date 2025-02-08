import logging
import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from scipy.sparse import coo_matrix
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def Umap_visualization(matrix, labels, pca) :
    umap_reducer = umap.UMAP(n_components=2, metric='euclidean', random_state=42)
    umap_results = umap_reducer.fit_transform(matrix)

    if pca == True:
        plt.subplot(1, 3, 1)
    else :
        plt.subplot(1, 2, 1)
    
    for cluster_label in set(labels):
        cluster_points = umap_results[labels == cluster_label]
        plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f"Cluster {cluster_label}")
    plt.legend()
    plt.title("UMAP Clustering Visualization")
    plt.xlabel("UMAP Dimension 1")
    plt.ylabel("UMAP Dimension 2")

def T_SNE_visualization(matrix, labels, pca) :
    tsne = TSNE(n_components=2, metric='euclidean', init="random", perplexity=30, random_state=42)
    tsne_results = tsne.fit_transform(matrix)

    if pca == True:
        plt.subplot(1, 3, 2)
    else :
        plt.subplot(1, 2, 2)
    for cluster_label in set(labels):
        if cluster_label == -1:
            cluster_points = tsne_results[labels == cluster_label]
            plt.scatter(cluster_points[:, 0], cluster_points[:, 1], c='gray', label="Noise", alpha=0.5)
        else:
            cluster_points = tsne_results[labels == cluster_label]
            plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f"Cluster {cluster_label}")
    plt.title("t-SNE Clustering Visualization")
    plt.xlabel("t-SNE Dimension 1")
    plt.ylabel("t-SNE Dimension 2")
    plt.legend()

def PCA_visualization(matrix, labels) :
    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(matrix)

    plt.subplot(1, 3, 3)
    unique_labels = set(labels)
    for label in unique_labels:
        cluster_points = reduced_data[labels == label]
        plt.scatter(
            cluster_points[:, 0],
            cluster_points[:, 1],
            label=f"Cluster {label}" if label != -1 else "Noise",
            s=50,
            alpha=0.7
        )
    plt.legend()
    plt.title("DBSCAN Clustering Visualization")
    plt.xlabel("PCA Dimension 1")
    plt.ylabel("PCA Dimension 2")

def visualization(matrix, labels, pca) :
    if pca :
        plt.figure(figsize=(20, 15))
        PCA_visualization(matrix, labels)
        Umap_visualization(matrix, labels, pca)
        T_SNE_visualization(matrix, labels, pca)
    else :
        plt.figure(figsize=(15, 15))
        Umap_visualization(matrix, labels, pca)
        T_SNE_visualization(matrix, labels, pca)
    plt.tight_layout()
    plt.show()


class Cluster :
    def __init__(self):
        self.matrix = None
        self.labels = None

    def clustering(self, file, sparse) :
        if not sparse :
            df = pd.read_csv(file, header=None, sep=' ')
            similarity_matrix = df.values

            similarity_matrix = np.clip(similarity_matrix, 0, 1)
            self.matrix = 1 - similarity_matrix

            dbscan = DBSCAN(eps=0.7, min_samples=3, metric='precomputed')
            self.labels = dbscan.fit_predict(self.matrix)

            print("Cluster Labels:", self.labels)
            return self
            
        else :
            data = []
            row = []
            col = []

            with open(file, 'r') as f:
                for line in f:
                    r, c, val = map(float, line.split())
                    distance_value = 1 - val
                    
                    row.append(int(r))
                    col.append(int(c))
                    data.append(distance_value)

            sparse_matrix = coo_matrix((data, (row, col))).tocsr()
            self.matrix = sparse_matrix

            dbscan = DBSCAN(eps=0.7, min_samples=2, metric='precomputed')
            self.labels = dbscan.fit_predict(self.matrix)

            print("Cluster Labels:", self.labels)
            return self

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    file_spect10_MH = "resultats_sim_cos/energy_10.0_precursor_M+H.txt"
    file_smil10_MH = "resultats_sim_smi/energy_10.0_precursor_M+H.txt"

    cluster = Cluster
    Cluster.clustering(cluster, file_smil10_MH, True)

    visualization(cluster.matrix, cluster.labels, False)