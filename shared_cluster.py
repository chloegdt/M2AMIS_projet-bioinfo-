import networkx as nx
from pyvis.network import Network

def load_clusters(file_path):
    
    clusters = {}
    with open(file_path, "r") as f:
        for cluster_id, line in enumerate(f):
            molecules = line.strip().split()
            clusters[cluster_id] = set(molecules)
    return clusters

def find_common_clusters(smiles_clusters, spectrum_clusters):
    
    common_clusters = {}

    for smiles_id, smiles_molecules in smiles_clusters.items():
        for spectrum_id, spectrum_molecules in spectrum_clusters.items():
            common_molecules = smiles_molecules & spectrum_molecules
            if common_molecules:
                common_clusters[(smiles_id, spectrum_id)] = list(common_molecules)

    return common_clusters

def plot_common_clusters(common_clusters, output_file="common_clusters.html"):
    
    G = nx.Graph()
    
    for (cluster_smi, cluster_spec), molecules in common_clusters.items():
        cluster_name = f"SMILES {cluster_smi} - Spectrum {cluster_spec}"
        G.add_node(cluster_name, color="red", size=20)
        
        for molecule in molecules:
            G.add_node(molecule, color="blue", size=10)
            G.add_edge(cluster_name, molecule, color="gray")

    if not G.nodes:
        print("⚠️ No common molecules found between clusters!")
        return
    
    net = Network(notebook=False, height="1000px", width="100%", bgcolor="black", font_color="white")
    
    for node, attrs in G.nodes(data=True):
        net.add_node(node, label=node, color=attrs.get("color", "gray"), size=attrs.get("size", 10))

    for edge in G.edges():
        net.add_edge(edge[0], edge[1], color="gray")
    net.write_html(output_file)

    print(f"Network saved as {output_file}. Open it manually in your browser.")


smiles_file = "resultats_clusters_dbscan/clusters_smiles_dbscan/energy_10.0_precursor_M+H.txt"
spectrum_file = "resultats_clusters_dbscan/clusters_spectres_dbscan/energy_10.0_precursor_M+H.txt"

file1 = ""
file2 = ""

smiles_clusters = load_clusters(smiles_file)
spectrum_clusters = load_clusters(spectrum_file)

common_clusters = find_common_clusters(smiles_clusters, spectrum_clusters)

plot_common_clusters(common_clusters)
