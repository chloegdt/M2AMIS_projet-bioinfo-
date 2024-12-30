import os
from matchms.importing import load_from_mgf
import matchms.filtering as ms_filters
import numpy as np
from matplotlib import pyplot as plt
from matchms import calculate_scores
from matchms.similarity import CosineGreedy

def metadata_processing(spectrum):
    spectrum = ms_filters.default_filters(spectrum)
    spectrum = ms_filters.repair_inchi_inchikey_smiles(spectrum)
    spectrum = ms_filters.derive_inchi_from_smiles(spectrum)
    spectrum = ms_filters.derive_smiles_from_inchi(spectrum)
    spectrum = ms_filters.derive_inchikey_from_inchi(spectrum)
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


path_data = "molecules_parsees/energy_102040_precursor_M-e.mgf"  # enter path to downloaded mgf file
file_mgf = os.path.join(path_data)
spectrums = list(load_from_mgf(file_mgf))

inchikeys = [s.get("inchikey") for s in spectrums]
found_inchikeys = np.sum([1 for x in inchikeys if x is not None])
print(f"Found {int(found_inchikeys)} inchikeys in metadata")

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
#plt.show()


# Je calcule la similarité des spectres avec la formule du cosinus
similarity_measure = CosineGreedy(tolerance=0.005)
scores = calculate_scores(spectrums, spectrums, similarity_measure, is_symmetric=True)
#print(scores.score_names)

# Je range les scores
scores_array = scores.scores.to_array()
scores_array[:5, :5]["CosineGreedy_score"]
#scores_array[:5, :5]["CosineGreedy_matches"]

#J'affiche le résultat. Les pics et les Smiles des molécules qui ont une valeur de Cosinus proche ou identique
best_matches = scores.scores_by_query(spectrums[5], name="CosineGreedy_score", sort=True)
print([x[1] for x in best_matches])
print([x[0].get("smiles") for x in best_matches])