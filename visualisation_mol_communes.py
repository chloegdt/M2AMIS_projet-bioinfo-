import networkx as nx
from pyvis.network import Network
from pyteomics import mgf

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

def molecules_pourcentage(file, cluster_spectre, cluster_smiles) :

    total_count = 0 
    with mgf.MGF(file) as spectres:
        for spectre in spectres:
            total_count += 1

    spectrum_molecules = set(mol for cluster in cluster_spectre.values() for mol in cluster)
    smiles_molecules = set(mol for cluster in cluster_smiles.values() for mol in cluster)

    shared_molecules = spectrum_molecules & smiles_molecules
    total_unique_molecules = spectrum_molecules | smiles_molecules

    percentage = (len(shared_molecules) / len(total_unique_molecules)) * 100 if total_unique_molecules else 0

    print(f"Nombre Total de molécules : ", total_count)
    print(f"Nombre de molécules dans les Clusters de Spectre : {len(spectrum_molecules)}")
    print(f"Nombre de molécules dans les Clusters de SMILES : {len(smiles_molecules)}")
    print(f"Molécules communes : {len(shared_molecules)}")
    print(f"Pourcentage de molécules communes : {percentage:.2f}%")
    

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

    print(f"Network saved as {output_file}. Ouvrez le fichier manuellement (firefox common_clusters).")

file = "energy_37.0_precursor_M+Na.mgf"
file1 = "cluster_smiles_37.0_precursor_M+Na.txt"
file2 = "cluster_spectre_37.0_precursor_M+Na.txt"

smiles_clusters = load_clusters(file1)
spectrum_clusters = load_clusters(file2)

common_clusters = find_common_clusters(smiles_clusters, spectrum_clusters)

plot_common_clusters(common_clusters)
molecules_pourcentage(file, spectrum_clusters, smiles_clusters)