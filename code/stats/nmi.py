import os
from sklearn.metrics import normalized_mutual_info_score

# à modifier selon les besoins
# BASE_DIR = os.getcwd()
BASE_DIR = "cluster_molecules"
HDBSCAN_DIR = os.path.join(BASE_DIR, "resultats_clusters_hdbscan")
MCL_DIR = os.path.join(BASE_DIR, "resultats_clusters_mcl")
DATA_TYPES = ["spectre", "fingerprint", "groupefonct"]
ENERGIES = ["10", "30", "50"]
OUTPUT_FILE = "res_nmi.txt"


def read_clusters(file_path):
    """
    Lit un fichier contenant les clusters et retourne un dictionnaire ID -> Cluster
    """
    cluster_dict = {}
    with open(file_path, 'r') as f:
        for cluster_id, line in enumerate(f):
            ids = line.strip().split("\t")
            for item in ids:
                cluster_dict[item] = cluster_id
    return cluster_dict


def compute_nmi(file1, file2):
    """
    Calcule la NMI entre deux fichiers de clustering
    """
    clusters1 = read_clusters(file1)
    clusters2 = read_clusters(file2)

    common_ids = set(clusters1.keys()).intersection(set(clusters2.keys()))
    if not common_ids:
        return None  # Pas d'IDs communs

    labels1 = [clusters1[id] for id in common_ids]
    labels2 = [clusters2[id] for id in common_ids]

    return normalized_mutual_info_score(labels1, labels2)

def main():
    """
    Parcourt les dossiers et compare tous les fichiers HDBSCAN vs MCL
    """
    base_dir = BASE_DIR
    hdbscan_dir = HDBSCAN_DIR
    mcl_dir = MCL_DIR

    data_types = DATA_TYPES
    energies = ENERGIES
    results = []


    for data_type in data_types:
    # Comparaisons HDBSCAN vs HDBSCAN au sein du même type de données
        for energy in energies:
            file_hdbscan1 = os.path.join(hdbscan_dir, data_type, f"energy_{energy}.txt")
            for data_type_2 in data_types:
                if data_type == data_type_2:  # Comparaison HDBSCAN vs HDBSCAN dans le même type
                    file_hdbscan2 = os.path.join(hdbscan_dir, data_type_2, f"energy_{energy}.txt")
                    relative_hdbscan_path1 = os.path.relpath(file_hdbscan1, base_dir)
                    relative_hdbscan_path2 = os.path.relpath(file_hdbscan2, base_dir)
                    if os.path.exists(file_hdbscan1) and os.path.exists(file_hdbscan2):
                        print(f"HDBSCAN{data_type} (Energy {energy}) vs HDBSCAN {data_type} (Energy {energy})")
                        nmi_score = compute_nmi(file_hdbscan1, file_hdbscan2)
                        if nmi_score is not None:
                            results.append(f"HDBSCAN {data_type} (Energy {energy}) vs HDBSCAN {data_type} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison HDBSCAN ({relative_hdbscan_path1}) vs ({relative_hdbscan_path2}) pour {data_type} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons MCL vs MCL au sein du même type de données
    for energy in energies:
        file_mcl1 = os.path.join(mcl_dir, data_type, f"energy_{energy}.txt")
        for data_type_2 in data_types:
            if data_type == data_type_2:  # Comparaison MCL vs MCL dans le même type
                file_mcl2 = os.path.join(mcl_dir, data_type_2, f"energy_{energy}.txt")
                relative_mcl_path1 = os.path.relpath(file_mcl1, base_dir)
                relative_mcl_path2 = os.path.relpath(file_mcl2, base_dir)
                if os.path.exists(file_mcl1) and os.path.exists(file_mcl2):
                    print(f"Comparaison MCL dans le même type : {relative_mcl_path1} vs {relative_mcl_path2}")
                    nmi_score = compute_nmi(file_mcl1, file_mcl2)
                    if nmi_score is not None:
                        results.append(f"MCL {data_type} (Energy {energy}) vs MCL {data_type} (Énergie {energy}): NMI = {nmi_score:.4f}")
                    else:
                        results.append(f"Comparaison MCL ({relative_mcl_path1}) vs ({relative_mcl_path2}) pour {data_type} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons HDBSCAN entre différents types de données
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison entre types différents
                    file_hdbscan1 = os.path.join(hdbscan_dir, data_type_1, f"energy_{energy}.txt")
                    file_hdbscan2 = os.path.join(hdbscan_dir, data_type_2, f"energy_{energy}.txt")
                    relative_hdbscan_path1 = os.path.relpath(file_hdbscan1, base_dir)
                    relative_hdbscan_path2 = os.path.relpath(file_hdbscan2, base_dir)
                    if os.path.exists(file_hdbscan1) and os.path.exists(file_hdbscan2):
                        print(f"Comparaison HDBSCAN entre différents types : {relative_hdbscan_path1} vs {relative_hdbscan_path2}")
                        nmi_score = compute_nmi(file_hdbscan1, file_hdbscan2)
                        if nmi_score is not None:
                            results.append(f"HDBSCAN {data_type_1} (Énergie {energy}) vs HDBSCAN {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison HDBSCAN ({relative_hdbscan_path1}) vs ({relative_hdbscan_path2}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons MCL entre différents types de données
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison entre types différents
                    file_mcl1 = os.path.join(mcl_dir, data_type_1, f"energy_{energy}.txt")
                    file_mcl2 = os.path.join(mcl_dir, data_type_2, f"energy_{energy}.txt")
                    relative_mcl_path1 = os.path.relpath(file_mcl1, base_dir)
                    relative_mcl_path2 = os.path.relpath(file_mcl2, base_dir)
                    if os.path.exists(file_mcl1) and os.path.exists(file_mcl2):
                        print(f"Comparaison MCL entre différents types : {relative_mcl_path1} vs {relative_mcl_path2}")
                        nmi_score = compute_nmi(file_mcl1, file_mcl2)
                        if nmi_score is not None:
                            results.append(f"MCL {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison MCL ({relative_mcl_path1}) vs ({relative_mcl_path2}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons HDBSCAN (même type) vs MCL (même type)
    for energy in energies:
        for data_type_1 in data_types:
            file_hdbscan = os.path.join(hdbscan_dir, data_type_1, f"energy_{energy}.txt")
            file_mcl = os.path.join(mcl_dir, data_type_1, f"energy_{energy}.txt")
            relative_hdbscan_path = os.path.relpath(file_hdbscan, base_dir)
            relative_mcl_path = os.path.relpath(file_mcl, base_dir)
            if os.path.exists(file_hdbscan) and os.path.exists(file_mcl):
                print(f"Comparaison HDBSCAN vs MCL pour le même type : {relative_hdbscan_path} vs {relative_mcl_path}")
                nmi_score = compute_nmi(file_hdbscan, file_mcl)
                if nmi_score is not None:
                    results.append(f"HDBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_1} (Énergie {energy}): NMI = {nmi_score:.4f}")
                else:
                    results.append(f"Comparaison HDBSCAN ({relative_hdbscan_path}) vs MCL ({relative_mcl_path}) pour {data_type_1} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons HDBSCAN (différents types) vs MCL (différents types)
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison HDBSCAN (différents types) vs MCL (différents types)
                    file_hdbscan = os.path.join(hdbscan_dir, data_type_1, f"energy_{energy}.txt")
                    file_mcl = os.path.join(mcl_dir, data_type_2, f"energy_{energy}.txt")
                    relative_hdbscan_path = os.path.relpath(file_hdbscan, base_dir)
                    relative_mcl_path = os.path.relpath(file_mcl, base_dir)
                    if os.path.exists(file_hdbscan) and os.path.exists(file_mcl):
                        print(f"HDBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy})")
                        nmi_score = compute_nmi(file_hdbscan, file_mcl)
                        if nmi_score is not None:
                            results.append(f"HDBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison HDBSCAN ({relative_hdbscan_path}) vs MCL ({relative_mcl_path}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")



    # Sauvegarde des résultats
    with open(OUTPUT_FILE, "w") as f:
        for line in results:
            f.write(line + "\n")

    print(f"Comparaisons terminées ! Résultats enregistrés dans {OUTPUT_FILE}.")

if __name__ == "__main__":
    main()

