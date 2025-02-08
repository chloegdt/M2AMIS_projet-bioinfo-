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

FILES = [
    #"energy_50.0_precursor_M+H.mgf",
    #"energy_10.0_precursor_M+H.mgf",
    #"energy_30.0_precursor_M+H.mgf",
    "energy_37.0_precursor_M+Na.mgf",
    "energy_25.0_precursor_M+Na.mgf"
]
MAPPING = {
    "visualisation",
    "similarite",
    "mcl",
    "parse",
    "fingerprint",
    "groups",
}
FILENAME1 = "Cluster.mgf"
FILENAME2 = "ALL_GNPS_cleaned.mgf"
DIRECTORY = "cluster_molecules/"


def get_parser():
    parser = argparse.ArgumentParser(
        description="Clustering of Spectrum and Smiles\nFichiers nécessaire :\nCluster.mgf et/ou ALL_GNPS_cleaned.mgf ",
        epilog="Commandes disponibles :\n"
                    "parse        - Vérifie et analyse les fichiers d'entrée.\n"
                    "similarite   - Calcule la similarité cosinus entre les spectres\n"
                    "               et la similarité fingerprints des smiles.\n"
                    "spectres     - Analyse les spectres de masse.\n"
                    "groups       - Regroupe les molécules en fonction de critères spécifiques.\n"
                    "functionnal  - Identifie les groupes fonctionnels présents.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "command",
        type=str,
        nargs="?",
        help="The command to execute (e.g., 'parse', 'similarite', etc.)."
    )
    parser.add_argument(
        "-p", "--path",
        type=pathlib.Path,
        help="Path to the input file or directory."
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

    elif command == "similarite" :
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
            cosinus.main_selected_files(DIRECTORY, FILES)

    elif command == "mcl" :
        inputdir = "cluster_molecules/resultats_cosinus_spectres"
        outputdir = "cluster_molecules/clusters_spectres_cosinus"
        markov_clustering_micans.clustering(inputdir, outputdir, "1.1")

    elif command == "visualisation" :
        print()

    elif command == "fingerprint":
        logging.info("Utilisation des fichiers par défaut.")
        logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
        check_chosen_files(parser)
        fingerprint.main(FILES)

    elif command == "groups":
        logging.info("Utilisation des fichiers par défaut.")
        logging.info("Pour utiliser d'autres fichiers : main.py {command} -p/--path PATH")
        check_chosen_files(parser)
        functionnal_group.main(FILES)

    else:
        parser.error(f"Commande inconnu: {command}. Utilisez 'help' pour voir les options disponibles.")



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    main()
