import matplotlib.pyplot as plt
from pyteomics import mgf
import numpy as np
import re

filename = "molecules_parsees/energy_102040_precursor_M-e.mgf"

# Je cherche le nombre de spectre dans le fichier
with mgf.MGF(filename) as spectres:
    spectre_list = list(spectres)
    num_spectre = len(spectre_list)

# Je fait une carte des couleurs basés sur le nombre de spectre
colormap = plt.colormaps.get_cmap('gist_rainbow').resampled(num_spectre)

plt.figure(figsize=(10, 6))

# Pour chaque IONS dans le fichier mgf je récupère la FORMULA et c'est valeurs m/z et leur intensité que j'affiche avec la couleur
for id, spectre in enumerate(spectre_list):
    color = colormap(id)
    formula = spectre['params'].get('formula')
    if formula == None or len(formula) <= 1:
        formula = "Unknown"
        print("HERE")
    

    mz = spectre['m/z array']
    intensity = spectre['intensity array']

    normalized_array = (intensity - intensity.min()) / (intensity.max() - intensity.min())

    plt.vlines(mz, ymin=0, ymax=normalized_array, color=color, linewidth=1, label=formula)


col_energy = re.findall(r'\d+', filename)[0]  # Get the first number

# Je fais le matplot pour afficher le résultat
plt.legend(title="Spectra", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(f"Mass Spectrum of collision energy: {col_energy}")
plt.xlabel("m/z")
plt.ylabel("Intensity")
plt.tight_layout()
plt.show()
