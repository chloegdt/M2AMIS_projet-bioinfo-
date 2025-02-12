import os
import multiprocessing
import logging
from scipy.sparse import coo_matrix
import numpy as np
import cosinus
from pathlib import Path
from matchms.importing import load_from_mgf
from matchms import calculate_scores
from matchms.filtering import normalize_intensities, select_by_intensity
from matchms.Spectrum import Spectrum
from matchms.similarity.BaseSimilarity import BaseSimilarity

class ManhattanSimilarity(BaseSimilarity):
    """
    Calcule la similarité Manhattan entre deux spectres.
    """
    score_datatype = np.dtype([('score', 'f4'), ('matched_peaks', 'i4')])

    def pair(self, reference: Spectrum, query: Spectrum):
        if reference is query:
            return np.array((1.0, len(np.union1d(reference.peaks.mz, query.peaks.mz))), dtype=self.score_datatype)
        
        ref_mz, ref_intensities = reference.peaks.mz, reference.peaks.intensities
        query_mz, query_intensities = query.peaks.mz, query.peaks.intensities
        common_mz = np.union1d(ref_mz, query_mz)
        ref_interp = np.interp(common_mz, ref_mz, ref_intensities, left=0, right=0)
        query_interp = np.interp(common_mz, query_mz, query_intensities, left=0, right=0)
        manhattan_distance = np.sum(np.abs(ref_interp - query_interp))
        similarity = 1 / (1 + manhattan_distance)
        matched_peaks = len(common_mz)
        return np.array((similarity, matched_peaks), dtype=self.score_datatype)

def save_scores(scores_data, output_file):
    """
    Sauvegarde les scores dans un fichier.
    """
    directory = os.path.dirname(output_file)
    if not os.path.exists(directory):
        os.makedirs(directory)
    rows = [item[0] for item in scores_data]
    cols = [item[1] for item in scores_data]
    scores = [item[2] for item in scores_data]
    sparse_matrix = coo_matrix((scores, (rows, cols)))
    with open(output_file, 'w') as f:
        for i, j, v in zip(sparse_matrix.row, sparse_matrix.col, sparse_matrix.data):
            f.write(f"{i+1} {j+1} {v}\n")

def compute_manhattan_scores(spectres):
    """
    Calcule la similarité Manhattan pour une liste de spectres.
    """
    similarity_measure = ManhattanSimilarity()
    scores = calculate_scores(spectres, spectres, similarity_measure, is_symmetric=False)
    scores_data = [(query.metadata.get("id") - 1, ref.metadata.get("id") - 1, score[0]) for query, ref, score in scores if query.metadata.get("id") <= ref.metadata.get("id")]
    return scores_data

def process_file(file_path, output_directory):
    """
    Traite un fichier MGF et calcule les scores Manhattan.
    """
    logging.info(f"Traitement du fichier : {file_path}")
    spectres = list(load_from_mgf(file_path))
    if not spectres:
        logging.error(f"Aucun spectre trouvé dans : {file_path}")
        return
    spectres = [cosinus.process(s) for s in spectres if s]
    spectres = cosinus.assign_ids_to_spectra(spectres)
    scores_data = compute_manhattan_scores(spectres)
    output_file = os.path.join(output_directory, f"{Path(file_path).stem}_manhattan_scores.txt")
    save_scores(scores_data, output_file)

def main(input_directory, output_directory):
    """
    Parcourt les fichiers d'un répertoire et applique la similarité Manhattan.
    """
    if not Path(input_directory).exists():
        logging.error(f"Répertoire introuvable : {input_directory}")
        return
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    files = [f for f in os.listdir(input_directory) if f.startswith("energy") and f.endswith(".mgf")]
    if not files:
        logging.error("Aucun fichier trouvé.")
        return
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        tasks = [(os.path.join(input_directory, file), output_directory) for file in files]
        pool.starmap(process_file, tasks)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    input_directory = "cluster_molecules_test"
    output_directory = "similarite_results"
    main(input_directory, output_directory)

