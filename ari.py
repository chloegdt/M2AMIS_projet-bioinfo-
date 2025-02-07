import os
import sys
from sklearn.metrics import adjusted_rand_score

def read_clusters(file_path):
    """
    Lit un fichier contenant les clusters et retourne un dictionnaire ID -> Cluster.
    """
    cluster_dict = {}
    with open(file_path, 'r') as f:
        for cluster_id, line in enumerate(f):
            ids = line.strip().split("\t")
            for item in ids:
                cluster_dict[item] = cluster_id
    return cluster_dict

def compute_ari(file1, file2):
    """
    Calcule l'ARI entre deux fichiers de clustering.
    """
    clusters1 = read_clusters(file1)
    clusters2 = read_clusters(file2)

    common_ids = set(clusters1.keys()).intersection(set(clusters2.keys()))
    if not common_ids:
        return None  # Aucun ID commun

    labels1 = [clusters1[id] for id in common_ids]
    labels2 = [clusters2[id] for id in common_ids]

    return adjusted_rand_score(labels1, labels2)

def main():
    """
    Compare tous les fichiers de deux dossiers et enregistre les ARI calculés.
    """
    if len(sys.argv) != 4:
        print("Utilisation : python ari.py <dossier1> <dossier2> <fichier_sortie>")
        sys.exit(1)

    dir1, dir2, output_file = sys.argv[1], sys.argv[2], sys.argv[3]

    # Vérifie que ce sont bien des dossiers
    if not (os.path.isdir(dir1) and os.path.isdir(dir2)):
        print("Les deux premiers arguments doivent être des dossiers valides.")
        sys.exit(1)

    # Récupère les fichiers .txt des deux dossiers
    files1 = {f for f in os.listdir(dir1) if f.endswith(".txt")}
    files2 = {f for f in os.listdir(dir2) if f.endswith(".txt")}

    # Intersection des fichiers présents dans les deux dossiers
    common_files = files1.intersection(files2)

    if not common_files:
        print("Aucun fichier commun trouvé entre les deux dossiers.")
        sys.exit(1)

    results = []

    for filename in common_files:
        file1 = os.path.join(dir1, filename)
        file2 = os.path.join(dir2, filename)

        print(f"Comparaison : {filename}")
        ari_score = compute_ari(file1, file2)

        if ari_score is not None:
            results.append(f"{filename} : ARI = {ari_score:.4f}")
        else:
            results.append(f"{filename} : Aucun ID commun")

    # Sauvegarde des résultats dans un fichier
    with open(output_file, "w") as f:
        for line in results:
            f.write(line + "\n")

    print(f"Comparaisons terminées ! Résultats enregistrés dans {output_file}.")

if __name__ == "__main__":
    main()

