import tkinter as tk
import matplotlib
matplotlib.use('TkAgg')  # Utiliser un backend TkAgg pour l'interactivité avec Tkinter

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Charger les clusters depuis le fichier texte
clusters = []
with open('clusters_smiles_tanimoto_I140/out.fg_energy_10.0_precursor_M+H.txt.I140', 'r') as file:
    for line in file:
        cluster = line.strip().split('\t')  # Diviser la ligne par tabulation
        clusters.append(cluster)

# Créer un graphe vide
G = nx.Graph()

# Ajouter les nœuds et les arêtes dans le graphe
for cluster in clusters:
    # Ajouter les nœuds dans le graphe
    G.add_nodes_from(cluster)
    
    # Ajouter des arêtes entre tous les nœuds du cluster pour les relier
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):
            G.add_edge(cluster[i], cluster[j])

# Créer une disposition des nœuds, en utilisant des positions adaptées aux clusters
# Utilisation d'une disposition de type "spring" (mais personnalisée pour les clusters)
pos = {}
x_offset = 0  # Décalage de position pour chaque cluster sur l'axe X

for cluster in clusters:
    # Positionner les nœuds d'un même cluster sur une zone spécifique
    y_offset = np.random.uniform(-1, 1)  # Aléatoire pour la position sur Y
    for i, node in enumerate(cluster):
        pos[node] = (x_offset + i * 0.1, y_offset)  # Disposition linéaire, mais un peu aléatoire
    x_offset += len(cluster) * 0.2  # Déplacer l'axe X pour chaque cluster

# Visualiser le graphe
plt.figure(figsize=(10, 8))
nx.draw(G, pos, with_labels=True, node_size=50, font_size=12, node_color="skyblue", edge_color="gray")
plt.show()
