import os
from matchms.importing import load_from_mgf
import matchms.filtering as ms_filters
import numpy as np
from matplotlib import pyplot as plt
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
import pandas as pd

def metadata_processing(spectrum):
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.repair_inchi_inchikey_smiles(spectrum)
    spectrum = ms_filters.harmonize_undefined_smiles(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchi(spectrum)
    spectrum = ms_filters.harmonize_undefined_inchikey(spectrum)
    spectrum = ms_filters.add_precursor_mz(spectrum)
    return spectrum

def peak_processing(spectrum):
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.normalize_intensities(spectrum)
    spectrum = ms_filters.select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = ms_filters.select_by_mz(spectrum, mz_from=10, mz_to=1000)
    return spectrum

def show_plot(spectrums) :
    numbers_of_peaks = [len(s.peaks.mz) for s in spectrums]

    # Filtrage des pics < 0.01
    spectrums = [metadata_processing(s) for s in spectrums]
    spectrums = [peak_processing(s) for s in spectrums]

    numbers_of_peaks2 = [len(s.peaks.mz) for s in spectrums]

    plt.figure(figsize=(6, 5), dpi=150)

    # Avant filtrage
    plt.hist(numbers_of_peaks, bins=20, alpha=0.5, color='blue', edgecolor='white', label='Before Filtering')

    # Après filtrage des pics < 0.01
    plt.hist(numbers_of_peaks2, bins=20, alpha=0.5, color='green', edgecolor='white', label='After Filtering')

    plt.title("Peaks per Spectrum (Before and After Filtering)")
    plt.xlabel("Number of Peaks in Spectrum")
    plt.ylabel("Number of Spectra")
    plt.legend()
    plt.tight_layout()
    plt.show()


def print_scores() :
    for (reference, query, score) in scores:
        print(f"Cosine score between {reference.get('formula')} and {query.get('formula')}" +
              f" is {score[0]:.2f} with {score[1]} matched peaks")


def save_to_csv(scores) :

    data = []
    for (reference, query, score) in scores:
        data.append({"reference": reference.metadata.get("formula", "unknown"),
                    "query": query.metadata.get("formula", "unknown"),
                    "score": score})

    df_scores = pd.DataFrame(data)

    output_csv_path = "similarity_scores.csv"
    df_scores.to_csv(output_csv_path, index=False)

    print(f"Similarity scores saved to {output_csv_path}")



####### MAIN ########
path_data = "molecules_parsees/energy_102040_precursor_M+Na.mgf"
file_mgf = os.path.join(path_data)
spectrums = list(load_from_mgf(file_mgf))

inchikeys = [s.get("inchikey") for s in spectrums]
found_inchikeys = np.sum([1 for x in inchikeys if x is not None])
print(f"Found {int(found_inchikeys)} inchikeys in metadata")


#show_plot(spectrums)

# Je calcule la similarité des spectres avec la formule du cosinus
similarity_measure = CosineGreedy(tolerance=0.005)
#scores = calculate_scores(spectrums, spectrums, similarity_measure, is_symmetric=True)
#print(scores.score_names)

scores = calculate_scores(spectrums, spectrums, CosineGreedy())

#print_scores()
save_to_csv(scores)