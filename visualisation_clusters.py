import logging
import networkx as nx
from pyvis.network import Network
import matplotlib.pyplot as plt
import numpy as np
import matplotlib


def read_clusters_from_file(file_path):
    """
        Permet de lire une matrice de cluster. Chaque ligne est un cluster
        et chaque nombre dans la ligne est l'id d'une molécule.

        @param Le chemin du fichier à lire.
        @return La liste des cluster.
    """
    clusters = {}
    
    with open(file_path, 'r') as f:
        for cluster_id, line in enumerate(f):
            molecule_ids = list(map(int, line.strip().split()))
            clusters[cluster_id] = molecule_ids
    
    return clusters




def read_cluster_matrix(file_path):
    """
        Permet de lire une matrice de cluster. Chaque ligne est un cluster
        et chaque nombre dans la ligne est l'id d'une molécule. Mais
        rajoute un nom à chaque cluster.

        @param Le chemin du fichier à lire.
        @return La liste des cluster.
    """
    clusters = {}

    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            values = line.strip().split()
            if values:
                cluster_name = f"Cluster {i+1}"
                clusters[cluster_name] = values

    return clusters




def plot_cluster_network(file_path):
    """
        Créer un plot des cluster. Mais la visualisation n'est pas très clair. 

        @param Le chemin du fichier à lire.
    """
    clusters = read_clusters_from_file(file_path)

    G = nx.Graph()

    num_clusters = len(clusters)
    cmap = matplotlib.colormaps.get_cmap('rainbow')
    color_norm = np.linspace(0, 1, num_clusters)
    cluster_colors = {cluster: cmap(color_norm[i]) for i, cluster in enumerate(clusters)}

    for cluster, molecules in clusters.items():
        cluster_name = f"Cluster {cluster}"
        G.add_node(cluster_name)
        
        for mol in molecules:
            molecule_name = f"Molecule {mol}"
            G.add_edge(cluster_name, molecule_name)

    
    node_colors = []
    for node in G.nodes():
        if "Cluster" in node:
            cluster_id = int(node.split()[1])
            node_colors.append(cluster_colors[cluster_id])
        else:
            node_colors.append("gray")
    
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    node_sizes = [2000 if "Cluster" in node else 500 for node in G.nodes()]
    nx.draw(G, pos, node_size=node_sizes, node_color=node_colors, edge_color="gray")

    plt.title("Graphe : Cluster-Molecules")
    plt.show()



def plot_interactive_network(file_path, output_file="network.html"):
    """
        Créer un fichier html qui est un network (un graphe). Dans ce graphe, chaque
        molécule est relié à son cluster.

        @param file_path : le chemin du fichier à lire.
        @param output_file="network.html" Le nom du fichier resultat.
    """

    clusters = read_cluster_matrix(file_path)

    G = nx.Graph()

    for cluster, molecules in clusters.items():
        G.add_node(cluster, color="red", size=20)
        for molecule in molecules:
            G.add_node(molecule, color="blue", size=10)
            G.add_edge(cluster, molecule, color="gray")

    if not G.nodes:
        print("No nodes found!")
        return

    net = Network(notebook=True, height="1000px", width="100%", bgcolor="white", font_color="white")
    
    for node, attrs in G.nodes(data=True):
        net.add_node(node, label=node, color=attrs.get("color", "gray"), size=attrs.get("size", 10))

    for edge in G.edges():
        net.add_edge(edge[0], edge[1], color="gray")

    try:
        net.show(output_file)
        #print(f"Saved as {output_file}")
        logging.info(f"Fichier de visualisation sauvergarder : {output_file}")
        logging.info(f"Pour ouvrir il faut faire : firefox {output_file}")
    except Exception as e:
        #print(f"Error network: {e}")
        logging.error(f"network: {e}")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    file = "test.txt"
    file_path = "clusters_smiles_dbscan/energy_10.0_precursor_M+H.txt"
    #plot_cluster_network(file_path)
    plot_interactive_network(file_path)
