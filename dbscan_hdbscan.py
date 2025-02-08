import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
import hdbscan

from scipy.sparse import coo_matrix
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn.utils.deprecation")
warnings.filterwarnings("ignore", category=UserWarning, module="umap")


FILE = "cluster_molecules/resultats_fingerprints/energy_37.0_precursor_M+Na.txt"


class DensityClustering:
    """
    Classe pour effectuer un clustering basé sur la densité (DBSCAN / HDBSCAN)
    à partir d'une matrice de similarité.
    """

    def __init__(self, file):
        """
        Initialise l'objet avec le chemin du fichier.

        Paramètres :
        ------------
        file : str
            Chemin vers le fichier de similarité entre spectres.
        """
        self.file = file
        self.matrix = None
        self.distance_matrix = None
        self.labels = None

    def load_data(self):
        """
        Charge les données de similarité à partir d'un fichier et construit
        une matrice creuse (sparse matrix).

        Format du fichier attendu :
        Chaque ligne contient :
            <id_spectre1> <id_spectre2> <similarité_cosinus>
        """
        rows, cols, data = [], [], []

        with open(self.file, 'r') as f:
            for line in f:
                spectre1, spectre2, similarity = line.split()
                spectre1, spectre2, similarity = int(spectre1) - 1, int(spectre2) - 1, float(similarity)
                rows.append(spectre1)
                cols.append(spectre2)
                data.append(similarity)

        self.matrix = coo_matrix((data, (rows, cols)))
        self.matrix = self.matrix.maximum(self.matrix.T)  # Garantit la symétrie
        self.distance_matrix = 1 - self.matrix.toarray()  # Conversion en matrice pleine

    def find_epsilon(self, min_samples=4):
        """
        Estime une valeur optimale d'epsilon pour DBSCAN en utilisant la courbe des distances.

        @param min_samples : int
            Nombre minimum de voisins pour considérer un point comme noyau.
        """
        neigh = NearestNeighbors(n_neighbors=min_samples, metric='precomputed')
        nbrs = neigh.fit(self.distance_matrix)
        distances, _ = nbrs.kneighbors(self.distance_matrix)
        distances = np.sort(distances[:, -1])  # Dernier voisin dans chaque ligne

        plt.plot(distances)
        plt.title('Valeur optimale d\'epsilon pour DBSCAN')
        plt.xlabel('Points triés')
        plt.ylabel('Distance au k-ième voisin')
        plt.show()

    def perform_clustering(self, clustering="HDBSCAN", eps=0.2, min_samples=4, metric="euclidean", reduction=True):
        """
        Applique le clustering DBSCAN ou HDBSCAN sur les données.

        @param clustering : str, par défaut "HDBSCAN", choix de l'algorithme ("DBSCAN" ou "HDBSCAN")
        @param eps : float, optionnel (DBSCAN uniquement) distance seuil pour regrouper des points
        @param in_samples : int, optionnel, par défaut 4, nombre minimal de voisins pour former un cluster
        @param reduction : boolean, optionnel, par défaut True, choix d'utiliser la réduction de dimension umap ou non
        @return np.ndarray : Labels des clusters attribués aux spectres.
        """
        if clustering == "HDBSCAN":
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples, metric=metric, alpha=1.0)
        elif clustering == "DBSCAN":
            clusterer = DBSCAN(eps=eps, min_samples=min_samples, metric=metric)
        else:
            print("ERREUR : ni DBSCAN ni HDBSCAN")
            return self.labels

        if reduction:
            reducer = umap.UMAP(n_components=2, metric='precomputed', random_state=42)
            umap_embedding = reducer.fit_transform(self.distance_matrix)
            self.labels = clusterer.fit_predict(umap_embedding)
        else:
            self.labels = clusterer.fit_predict(self.distance_matrix)

        # print("Labels de clusters :")
        # print(np.unique(self.labels, return_counts=True))
        n_clusters = len(set(self.labels)) - (1 if -1 in self.labels else 0)
        logging.info(f"Nombre de clusters détectés : {n_clusters}")
        return self.labels

    def visualize_clusters(self, clustering="HDBSCAN"):
        """
        Affiche une visualisation des clusters en 2D après réduction de dimension avec UMAP.

        @param clustering : str, optionnel, nom du clustering utilisé (pour le titre du graphique)
        """
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
        """
        Enregistre les clusters obtenus dans un fichier

        Format du fichier :
        Chaque ligne représente un cluster :
            <id_spectre1> <id_spectre2> ...

        @param output_file : str, chemin du fichier de sortie
        """
        clusters_dict = {}
        for idx, label in enumerate(self.labels):
            if label != -1:
                clusters_dict.setdefault(label, []).append(idx + 1)

        with open(output_file, 'w') as f:
            for cluster_id, ids in clusters_dict.items():
                f.write("\t".join(map(str, ids)) + "\n")



def clustering_hdbscan(files, input_directory, output_directory, reduction=True):
    """
    Applique le clustering HDBSCAN sur les fichiers d'entrée et sauvegarde en fichier texte.
    
    @param files: Liste contenant les nom des fichiers à traiter.
    @param input_directory: chemin du dossier des fichers d'entrées.
    @param output_directory: chemin du dossier de sauvegarde.
    @param reduction: booleen optionnel, choix d'utiliser la réduction de dimension umap ou non.
    """
    if reduction:
        min_samples = 4
        metric = "euclidean"
    else:
        min_samples = 30
        metric = "precomputed"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory, exist_ok=True)
        logging.info(f"Répertoire {output_directory} créé.")

    for file in files:
        file_path = os.path.join(input_directory, f"{os.path.splitext(os.path.basename(file))[0]}.txt")
        if not os.path.isfile(file_path):
            logging.warning(f"Fichier {file_path} introuvable.")
            continue

        logging.info(f"Clustering HDBSCAN du fichier {file}.")
        clustering = DensityClustering(file_path)
        clustering.load_data()
        labels = clustering.perform_clustering("HDBSCAN", min_samples=min_samples, metric=metric, reduction=reduction)
        clustering.save_clusters_to_file(os.path.join(output_directory, os.path.basename(file_path)))
    logging.info(f"Clusters sauvegardés dans le dossier : {output_directory}\n")



if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    # file = FILE

    # clustering = DensityClustering(file)
    # clustering.load_data()
    # labels = clustering.perform_clustering("DBSCAN", eps=0.2, min_samples=4)
    # clustering.save_clusters_to_file(OUTPUT)
    # clustering.find_epsilon(4)
    # print(f"Bruit : {sum(1 for l in labels if l == -1)} points")
    # clustering.visualize_clusters("DBSCAN")
    clustering_hdbscan([FILE], "cluster_molecules/resultats_fingerprints/", "cluster_molecules/clusters_fingerprints/", True)
