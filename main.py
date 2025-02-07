import os
import sys
import logging
import argparse
import pathlib

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
    "energy_25.0_precursor_M+Na.mgf",
]

MAPPING = {
    #"fingerprints": fingerprint.main,
    #"fg": fingerprint.main,
    #"functionnal": functionnal_group.main,
    #"groups": functionnal_group.main,
    #"dbscan": dbscan_hdbscan.main,
    "similarite",
    #"markov": markov_clustering_micans.main,
    #"mcl": markov_clustering_micans.main,
    "parse"
}

def get_parser():

    parser = argparse.ArgumentParser(
        description="Clustering of Spectrum and Smiles",
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
    directory = "cluster_molecules/"
    for file in FILES:
        if not os.path.exists(os.path.join(directory, file)):
            print("Les fichiers nécessaires ne sont pas présents dans", directory, "\nExécution de parse...")
            parse_data_final.main()
            break


def main():
    parser = get_parser()
    args = parser.parse_args()

    if args.command is None or args.command.lower() == "help":
        parser.print_help()
        sys.exit(0)
    
    command = args.command.lower()

    if command == "parse" :
        check_files()

    elif command == "similarite" :
        if args.path:
            print(f"Path provided: {args.path.resolve()}")
            if not args.path.exists():
                print("Error: Path does not exist.")
                sys.exit(1)
        else:
            print("No path provided. Using preconfigured files.")
            print("If you want to use a directory : \nrequired : main.py similarite -p, --path PATH")
            directory_path = "cluster_molecules_chosen"
            cosinus.main_selected_files(directory_path, FILES)

        
    else:
        parser.error(f"Code inconnu: {command}. Utilisez 'help' pour voir les options disponibles.")


if __name__ == "__main__":
    #logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    main()