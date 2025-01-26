import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_similarity_data(file_path):
    """
    Charge les scores de similarité depuis un fichier CSV.
    """
    try:
        data = pd.read_csv(file_path)
       # Nettoyer la colonne "score" en extrayant les valeurs
        data[['similarity', 'matched_peaks']] = data['score'].str.strip("[]").str.split(", ", expand=True)
        
        # Convertir les nouvelles colonnes en float/int
        data['similarity'] = data['similarity'].astype(float)
        data['matched_peaks'] = data['matched_peaks'].astype(int)
        return data
    except FileNotFoundError:
        print(f"Erreur : le fichier {file_path} est introuvable.")
        return None
    except Exception as e:
        print(f"Erreur lors du nettoyage des données : {e}")
        return None

def plot_similarity_heatmap(data, output_path="heatmap.png"):
    """
    Crée une matrice de similarité sous forme de heatmap.
    """
    # Identifier les spectres uniques
    spectrums = sorted(set(data['reference']).union(set(data['query'])))
    
    # Construire une matrice de similarité
    similarity_matrix = pd.DataFrame(
        np.zeros((len(spectrums), len(spectrums))),
        index=spectrums,
        columns=spectrums
    )
    
    # Remplir la matrice avec les scores de similarité
    for _, row in data.iterrows():
        similarity_matrix.at[row['reference'], row['query']] = row['score']
        similarity_matrix.at[row['query'], row['reference']] = row['score']  # Symétrie

    # Dessiner la heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_matrix, cmap="YlGnBu", square=True, cbar_kws={'label': 'Similarité Manhattan'})
    plt.title("Matrice de Similarité Manhattan")
    plt.xlabel("Spectres")
    plt.ylabel("Spectres")
    plt.tight_layout()
    
    # Sauvegarder la heatmap
    plt.savefig(output_path)
    print(f"Heatmap sauvegardée sous : {output_path}")
    plt.show()

#### MAIN ####
if __name__ == "__main__":
    # Chemin vers le fichier CSV généré par votre programme principal
    input_csv_path = "similarite_results/energy_0.0_precursor_M+H_manhattan_scores.csv"  # Remplacez par le chemin réel
    output_heatmap_path = "manhattan_heatmap.png"

    # Charger les données
    similarity_data = load_similarity_data(input_csv_path)

    if similarity_data is not None:
        # Créer la heatmap
        plot_similarity_heatmap(similarity_data, output_path=output_heatmap_path)
