import os
import sys
import logging
import argparse

import parse_data_final
import cosinus
import fingerprint
import functionnal_group
import dbscan_hdbscan
import markov_clustering_micans

FILES = [
    "energy_50.0_precursor_M+H.mgf",
    "energy_10.0_precursor_M+H.mgf",
    "energy_30.0_precursor_M+H.mgf",
    "energy_25.0_precursor_M+Na.mgf",
]

MAPPING = {
    #"fingerprints": fingerprint.main,
    #"fg": fingerprint.main,
    #"functionnal": functionnal_group.main,
    #"groups": functionnal_group.main,
    #"dbscan": dbscan_hdbscan.main,
    #"cosinus": cosinus.main,
    #"spectres": cosinus.main,
    #"markov": markov_clustering_micans.main,
    #"mcl": markov_clustering_micans.main,
    "parse": parse_data_final.main,
}

def get_parser():
    description = """\
Code du module à exécuter. Utilisez 'help' pour afficher la liste des codes disponibles.

Commandes disponibles :
  parse        - Vérifie et analyse les fichiers d'entrée.
  cosinus      - Calcule la similarité cosinus entre les spectres.
  spectres     - Analyse les spectres de masse.
  groups       - Regroupe les molécules en fonction de critères spécifiques.
  functionnal  - Identifie les groupes fonctionnels présents.
  fingerprints - Génère des empreintes moléculaires.
  fg           - Alias pour 'fingerprints', génère des empreintes moléculaires.
"""

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter  # Keeps the formatting intact
    )

    parser.add_argument(
        "code",
        type=str,
        nargs="?",
        help="Le module à exécuter (voir la liste ci-dessous)."
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

    if args.code is None or args.code.lower() == "help":
        parser.print_help()
        sys.exit(0)
    
    code = args.code.lower()

    if code == "parse" :
        check_files()
    elif code in ["cosinus", "spectres", "groups", "functionnal", "fingerprints", "fg"]:
        print()
    else:
        print("je suis aussi ici")
        parser.error(f"Code inconnu: {code}. Utilisez 'help' pour voir les options disponibles.")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    main()