import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats
from collections import Counter


FILE1 = os.path.join("cluster_molecules", "resultats_clusters_hdbscan","spectre","energy_30.txt")
FILE2 = os.path.join("cluster_molecules", "resultats_clusters_hdbscan","fingerprint","energy_30.txt")
OUTPUT_FILE = "comparison_results_spectre1_smiles2_30_pourplot.txt"
HISTOGRAM1 =  "histogram_spectre30.png"
TITRE1 = "Distribution des tailles des clusters de nos spectres"
HISTOGRAM2 = "histogram_smiles30.png"
TITRE2 = "Distribution des tailles des clusters de nos FM"


def analyze_clusters(file_path):
    """
    Récupère les clusters depuis le fichier de résultat.
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    clusters = [set(line.strip().split('\t')) for line in lines if line.strip()]
    return clusters

def jaccard_index(c1, c2):
    return len(c1 & c2) / len(c1 | c2) if c1 | c2 else 0

def compute_overlap(clusters):
    overlaps = [len(c1 & c2) / len(c1) for c1 in clusters for c2 in clusters if c1 != c2]
    return np.mean(overlaps) if overlaps else 0

def save_results(stats, output_file):
    """
    Enregistre les résultats dans un fichier texte
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")


def plot_histogram(sizes, title, output_file):
    plt.figure()
    plt.hist(sizes, bins=50, alpha=0.7, color='blue', edgecolor='black')
    plt.xlim(0, 85)
    plt.xlabel("Taille des clusters")
    plt.ylabel("Nombre de clusters")
    plt.title(title)
    plt.savefig(output_file)
    plt.close()


def compare_clusters(file1, file2):
    """
    Compare deux fichiers de clusters et génère des statistiques.
    """
    clusters1 = analyze_clusters(file1)
    clusters2 = analyze_clusters(file2)
    
    num_clusters1 = len(clusters1)
    num_clusters2 = len(clusters2)
    
    common_clusters = set(map(frozenset, clusters1)) & set(map(frozenset, clusters2))
    unique_clusters1 = set(map(frozenset, clusters1)) - common_clusters
    unique_clusters2 = set(map(frozenset, clusters2)) - common_clusters
    
    all_ids1 = set().union(*clusters1)
    all_ids2 = set().union(*clusters2)
    common_ids = all_ids1 & all_ids2
    unique_ids1 = all_ids1 - all_ids2
    unique_ids2 = all_ids2 - all_ids1
    
    subsets1_in_2 = [(len(c1), len(c2), len(c1 & c2)) for c1 in clusters1 for c2 in clusters2 if c1.issubset(c2)]
    subsets2_in_1 = [(len(c2), len(c1), len(c2 & c1)) for c1 in clusters1 for c2 in clusters2 if c2.issubset(c1)]
    
    sizes1 = [len(cluster) for cluster in clusters1]
    sizes2 = [len(cluster) for cluster in clusters2]
    
    overlap_rate1 = compute_overlap(clusters1)
    overlap_rate2 = compute_overlap(clusters2)
    
    id_cluster_count1 = {id_: sum(1 for cluster in clusters1 if id_ in cluster) for id_ in all_ids1}
    id_cluster_count2 = {id_: sum(1 for cluster in clusters2 if id_ in cluster) for id_ in all_ids2}
    avg_id_cluster_count1 = np.mean(list(id_cluster_count1.values())) if id_cluster_count1 else 0
    avg_id_cluster_count2 = np.mean(list(id_cluster_count2.values())) if id_cluster_count2 else 0
    
    jaccard_distances = [jaccard_index(c1, c2) for c1 in clusters1 for c2 in clusters2 if c1 != c2]
    avg_jaccard_distance = np.mean(jaccard_distances) if jaccard_distances else 0
    
    stats = {
        "num_clusters1": num_clusters1,
        "num_clusters2": num_clusters2,
        "common_clusters": len(common_clusters),
        "unique_clusters1": len(unique_clusters1),
        "unique_clusters2": len(unique_clusters2),
        "common_ids": len(common_ids),
        "unique_ids1": len(unique_ids1),
        "unique_ids2": len(unique_ids2),
        "subsets1_in_2": subsets1_in_2,
        "subsets2_in_1": subsets2_in_1,
        "avg_size1": np.mean(sizes1) if sizes1 else 0,
        "avg_size2": np.mean(sizes2) if sizes2 else 0,
        "median_size1": np.median(sizes1) if sizes1 else 0,
        "median_size2": np.median(sizes2) if sizes2 else 0,
        "std_size1": np.std(sizes1) if sizes1 else 0,
        "std_size2": np.std(sizes2) if sizes2 else 0,
        "max_size1": max(sizes1, default=0),
        "min_size1": min(sizes1, default=0),
        "max_size2": max(sizes2, default=0),
        "min_size2": min(sizes2, default=0),
        "size_distribution1": dict(Counter(sizes1)),
        "size_distribution2": dict(Counter(sizes2)),
        "overlap_rate1": overlap_rate1,
        "overlap_rate2": overlap_rate2,
        "avg_id_cluster_count1": avg_id_cluster_count1,
        "avg_id_cluster_count2": avg_id_cluster_count2,
        "avg_jaccard_distance": avg_jaccard_distance
    }
    
    return stats, sizes1, sizes2

if __name__ == "__main__":
    file1 = FILE1
    file2 = FILE2
    output_file = OUTPUT_FILE
    
    stats, sizes1, sizes2 = compare_clusters(file1, file2)
    plot_histogram(sizes1, TITRE1, HISTOGRAM1)
    plot_histogram(sizes2, TITRE2, HISTOGRAM2)
    save_results(stats, output_file)
    print(f"Les résultats ont été enregistrés dans {output_file}")

