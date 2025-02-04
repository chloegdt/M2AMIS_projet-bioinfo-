import os
from sklearn.metrics import normalized_mutual_info_score

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
    Parcourt les dossiers et compare tous les fichiers DBSCAN vs MCL
    """
    base_dir = os.getcwd()
    dbscan_dir = os.path.join(base_dir, "resultats_clusters_dbscan")
    mcl_dir = os.path.join(base_dir, "resultats_clusters_mcl")

    data_types = ["spectres", "smiles", "descriptors"]
    energies = ["10", "30", "50"]

    results = []


    for data_type in data_types:
    # Comparaisons DBSCAN vs DBSCAN au sein du même type de données
        for energy in energies:
            file_dbscan1 = os.path.join(dbscan_dir, data_type, f"{data_type}_{energy}.txt")
            for data_type_2 in data_types:
                if data_type == data_type_2:  # Comparaison DBSCAN vs DBSCAN dans le même type
                    file_dbscan2 = os.path.join(dbscan_dir, data_type_2, f"{data_type_2}_{energy}.txt")
                    relative_dbscan_path1 = os.path.relpath(file_dbscan1, base_dir)
                    relative_dbscan_path2 = os.path.relpath(file_dbscan2, base_dir)
                    if os.path.exists(file_dbscan1) and os.path.exists(file_dbscan2):
                        print(f"DBSCAN{data_type} (Energy {energy}) vs DBSCAN {data_type} (Energy {energy})")
                        nmi_score = compute_nmi(file_dbscan1, file_dbscan2)
                        if nmi_score is not None:
                            results.append(f"DBSCAN {data_type} (Energy {energy}) vs DBSCAN {data_type} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison DBSCAN ({relative_dbscan_path1}) vs ({relative_dbscan_path2}) pour {data_type} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons MCL vs MCL au sein du même type de données
    for energy in energies:
        file_mcl1 = os.path.join(mcl_dir, data_type, f"{data_type}_{energy}.txt")
        for data_type_2 in data_types:
            if data_type == data_type_2:  # Comparaison MCL vs MCL dans le même type
                file_mcl2 = os.path.join(mcl_dir, data_type_2, f"{data_type_2}_{energy}.txt")
                relative_mcl_path1 = os.path.relpath(file_mcl1, base_dir)
                relative_mcl_path2 = os.path.relpath(file_mcl2, base_dir)
                if os.path.exists(file_mcl1) and os.path.exists(file_mcl2):
                    print(f"Comparaison MCL dans le même type : {relative_mcl_path1} vs {relative_mcl_path2}")
                    nmi_score = compute_nmi(file_mcl1, file_mcl2)
                    if nmi_score is not None:
                        results.append(f"MCL {data_type} (Energy {energy}) vs MCL {data_type} (Énergie {energy}): NMI = {nmi_score:.4f}")
                    else:
                        results.append(f"Comparaison MCL ({relative_mcl_path1}) vs ({relative_mcl_path2}) pour {data_type} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons DBSCAN entre différents types de données
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison entre types différents
                    file_dbscan1 = os.path.join(dbscan_dir, data_type_1, f"{data_type_1}_{energy}.txt")
                    file_dbscan2 = os.path.join(dbscan_dir, data_type_2, f"{data_type_2}_{energy}.txt")
                    relative_dbscan_path1 = os.path.relpath(file_dbscan1, base_dir)
                    relative_dbscan_path2 = os.path.relpath(file_dbscan2, base_dir)
                    if os.path.exists(file_dbscan1) and os.path.exists(file_dbscan2):
                        print(f"Comparaison DBSCAN entre différents types : {relative_dbscan_path1} vs {relative_dbscan_path2}")
                        nmi_score = compute_nmi(file_dbscan1, file_dbscan2)
                        if nmi_score is not None:
                            results.append(f"DBSCAN {data_type_1} (Énergie {energy}) vs DBSCAN {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison DBSCAN ({relative_dbscan_path1}) vs ({relative_dbscan_path2}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons MCL entre différents types de données
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison entre types différents
                    file_mcl1 = os.path.join(mcl_dir, data_type_1, f"{data_type_1}_{energy}.txt")
                    file_mcl2 = os.path.join(mcl_dir, data_type_2, f"{data_type_2}_{energy}.txt")
                    relative_mcl_path1 = os.path.relpath(file_mcl1, base_dir)
                    relative_mcl_path2 = os.path.relpath(file_mcl2, base_dir)
                    if os.path.exists(file_mcl1) and os.path.exists(file_mcl2):
                        print(f"Comparaison MCL entre différents types : {relative_mcl_path1} vs {relative_mcl_path2}")
                        nmi_score = compute_nmi(file_mcl1, file_mcl2)
                        if nmi_score is not None:
                            results.append(f"MCL {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison MCL ({relative_mcl_path1}) vs ({relative_mcl_path2}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons DBSCAN (même type) vs MCL (même type)
    for energy in energies:
        for data_type_1 in data_types:
            file_dbscan = os.path.join(dbscan_dir, data_type_1, f"{data_type_1}_{energy}.txt")
            file_mcl = os.path.join(mcl_dir, data_type_1, f"{data_type_1}_{energy}.txt")
            relative_dbscan_path = os.path.relpath(file_dbscan, base_dir)
            relative_mcl_path = os.path.relpath(file_mcl, base_dir)
            if os.path.exists(file_dbscan) and os.path.exists(file_mcl):
                print(f"Comparaison DBSCAN vs MCL pour le même type : {relative_dbscan_path} vs {relative_mcl_path}")
                nmi_score = compute_nmi(file_dbscan, file_mcl)
                if nmi_score is not None:
                    results.append(f"DBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_1} (Énergie {energy}): NMI = {nmi_score:.4f}")
                else:
                    results.append(f"Comparaison DBSCAN ({relative_dbscan_path}) vs MCL ({relative_mcl_path}) pour {data_type_1} (Énergie {energy}): AUCUN ID COMMUN")

    # Comparaisons DBSCAN (différents types) vs MCL (différents types)
    for energy in energies:
        for data_type_1 in data_types:
            for data_type_2 in data_types:
                if data_type_1 != data_type_2:  # Comparaison DBSCAN (différents types) vs MCL (différents types)
                    file_dbscan = os.path.join(dbscan_dir, data_type_1, f"{data_type_1}_{energy}.txt")
                    file_mcl = os.path.join(mcl_dir, data_type_2, f"{data_type_2}_{energy}.txt")
                    relative_dbscan_path = os.path.relpath(file_dbscan, base_dir)
                    relative_mcl_path = os.path.relpath(file_mcl, base_dir)
                    if os.path.exists(file_dbscan) and os.path.exists(file_mcl):
                        print(f"DBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy})")
                        nmi_score = compute_nmi(file_dbscan, file_mcl)
                        if nmi_score is not None:
                            results.append(f"DBSCAN {data_type_1} (Énergie {energy}) vs MCL {data_type_2} (Énergie {energy}): NMI = {nmi_score:.4f}")
                        else:
                            results.append(f"Comparaison DBSCAN ({relative_dbscan_path}) vs MCL ({relative_mcl_path}) pour {data_type_1} vs {data_type_2} (Énergie {energy}): AUCUN ID COMMUN")



    # Sauvegarde des résultats
    with open("resultats_nmi.txt", "w") as f:
        for line in results:
            f.write(line + "\n")

    print("Comparaisons terminées ! Résultats enregistrés dans 'resultats_nmi.txt'.")

if __name__ == "__main__":
    main()

