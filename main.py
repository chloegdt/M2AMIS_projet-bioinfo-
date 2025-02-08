import os
import sys
import logging
import argparse
import pathlib
from pathlib import Path

import parse_data_final
import cosinus
import fingerprint
import functionnal_group
import dbscan_hdbscan
import markov_clustering_micans
from visualisation_clusters import plot_interactive_network
from visualisation_mol_communes import find_common_clusters, load_clusters, net_common_clusters

FILES = [
    # "energy_50.0_precursor_M+H.mgf",
    # "energy_10.0_precursor_M+H.mgf",
    # "energy_30.0_precursor_M+H.mgf",
    "energy_37.0_precursor_M+Na.mgf",
    "energy_25.0_precursor_M+Na.mgf"
]
MAPPING = {
    "visualisation",
    "cosinus",
    "mcl",
    "parse",
    "fingerprint",
    "groupes",
    "hdbscan",
}
FILENAME1 = "Cluster.mgf"
FILENAME2 = "ALL_GNPS_cleaned.mgf"
DIRECTORY = "cluster_molecules/"


def get_parser():
    parser = argparse.ArgumentParser(
        description="Clustering of Spectrum and Smiles\nFichiers nécessaire :\nCluster.mgf et/ou ALL_GNPS_cleaned.mgf ",
        epilog="Commandes disponibles :\n"
                    "parse          - Vérifie et analyse les fichiers d'entrée.\n"
                    "cosinus        - Calcule la similarité cosinus entre les spectres\n"
                    "                 et la similarité fingerprints des smiles.\n"
                    "spectres       - Analyse les spectres de masse.\n"
                    "groupes        - Regroupe les molécules en fonction de critères spécifiques.\n"
                    "functionnal    - Identifie les groupes fonctionnels présents.\n"
                    "hdbscan        - Applique le clustering HDBSCAN\n"
                    "intersection   - Permet de visualiser les clusters commun entre deux fichiers de cluster\n"
                    "cluster        - Permet de visualiser les cluster d'un fichier.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "command",
        type=str,
        nargs="?",
        help="Les commandes à executé(ex : 'parse', 'cosinus', etc.)."
    )
    parser.add_argument(
        "-p", "--path",
        type=pathlib.Path,
        help="Chemin vers un fichier ou repertoire.\nNecessaire pour cluster et intersection."
    )
    parser.add_argument(
        "-s",
        type=pathlib.Path,
        help="Chemin vers un fichier.\nNecessaire pour intersection."
    )
    return parser


def check_files():
    for file in FILES:
        if not os.path.exists(os.path.join(DIRECTORY, file)):
            logging.warning(f"Les fichiers nécessaires ne sont pas présents dans {DIRECTORY} \nExécution de parse...")
            parse_data_final.main()
            break


def check_chosen_files(parser):
    count = 0
    for file in FILES :
        file = "cluster_molecules/" + file 
        if os.path.isfile(file) :
            count += 1

    if count == 0:
        print("Il faut que au moins 1 des 5 fichiers soit présent après le parsing.")
        parser.print_help()
        sys.exit(0)

def check_path(args):
    if args.path:
        if not args.path.exists():
            logging.error("Le fichier donné est introuvable.")
            return False
        print(f"File provided: {args.path.resolve()}")
        file_path = args.path
        print(file_path)
        return True
    else:
        logging.error("Il faut un chemin vers un fichier à visualiser.")
        return False

def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.command is None or args.command.lower() == "help":
        parser.print_help()
        sys.exit(0)
    if not os.path.isfile(FILENAME1) and not os.path.isfile(FILENAME2):
        logging.error(f"Aucun fichier .mgf trouvé ({FILENAME1} et/ou {FILENAME2}).")
        parser.print_help()
        sys.exit(0)

    command = args.command.lower()
    if command == "parse" :
        check_files()

    elif command == "cosinus" :
        if args.path:
            print(f"Path provided: {args.path.resolve()}")
            directory_path = args.path
            print(directory_path)
            cosinus.main(directory_path)
            if not args.path.exists():
                logging.error("Le chemin entré est introuvable.")
                sys.exit(1)
        else:
            logging.info("Utilisation des fichiers par défaut.")
            logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
            check_chosen_files(parser)
            cosinus.main_selected_files(FILES, DIRECTORY)

    elif command == "mcl" :
        inputdir = "cluster_molecules/resultats_cosinus_spectres"
        outputdir = "cluster_molecules/clusters_spectres_cosinus"
        markov_clustering_micans.clustering(inputdir, outputdir, "2.0")

    elif command == "intersection":
        if check_path(args) :
            smiles_clusters = load_clusters(args.path)
            spectrum_clusters = load_clusters(args.path)

            common_clusters = find_common_clusters(smiles_clusters, spectrum_clusters)

            net_common_clusters(common_clusters)
        else :
            logging.error("Fichier demandé : .txt de deux clusters")

    elif command == "cluster":
        if check_path(args) :
            plot_interactive_network(args.path)
        else :
            logging.error("Fichier demandé : .txt de cluster")
            logging.error("Exemple : cluster_molecules/HDBSCAN_cosinus_spectre/energy_37.0_precursor_M+Na.txt")

    elif command == "fingerprint":
        logging.info("Utilisation des fichiers par défaut.")
        logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
        check_chosen_files(parser)
        fingerprint.main(FILES)

    elif command == "groupes":
        logging.info("Utilisation des fichiers par défaut.")
        logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
        check_chosen_files(parser)
        functionnal_group.main(FILES)

    elif command == "hdbscan":
        logging.info("Utilisation des fichiers par défaut.")
        logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
        check_chosen_files(parser)
        logging.info("HDBSCAN sur spectres.")
        dbscan_hdbscan.clustering_hdbscan(FILES, "cluster_molecules/resultats_cosinus_spectres/", "cluster_molecules/HDBSCAN_cosinus_spectre/", True)
        logging.info("HDBSCAN sur fingerprints.")
        dbscan_hdbscan.clustering_hdbscan(FILES, "cluster_molecules/resultats_fingerprints/", "cluster_molecules/HDBSCAN_fingerprints_smiles/", True)
        logging.info("HDBSCAN sur groupes fonctionnels.")
        dbscan_hdbscan.clustering_hdbscan(FILES, "cluster_molecules/resultats_groupes-fonc/", "cluster_molecules/HDBSCAN_groupes_smiles/", False)

    else:
        parser.error(f"Commande inconnu: {command}. Utilisez 'help' pour voir les options disponibles.")



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    main()
