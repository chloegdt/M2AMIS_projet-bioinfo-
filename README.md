# Module Projet - M2 AMIS - Février 2025 
#### Clémence DUMOULIN, Sébastien FERNANDEZ, Chloé GODET,
#### Raphaël GUERIN Jeffried NOULEHO NOUADJE
________________________

## Description

Ce projet, réalisé dans le cadre du Master 2 AMIS, vise à étudier la similarité entre deux représentations de molécules : les spectres de masses et les SMILES.  
Le code est principalement écrit en Python et utilise diverses bibliothèques pour le traitement et l'analyse des données.


## Structure du projet

* analyse_clusters_groupesFonctionnels.py : Analyse les clusters de molécules selon leurs groupes fonctionnels dominant.
* cosinus.py : Traite les spectres des molécules puis calcule la similarité cosinus.
* dbscan_hdbscan.py : Applique les algorithmes de clustering DBSCAN et HDBSCAN sur des fichiers de similarité.
* fingerprint.py : Génère des fingerprints de Morgan et calcule la similarité de Tanimoto.
* functionnal_group.py : Identifie la présence de groupes fonctionnels et calcule la similarité de Tanimoto.
* main.py : Gère l'éxecution des différents programmes.
* markov_clustering_micans.py : Réalise le clustering des molécules en utilisant l'algorithme de Markov Clustering (MCL).
* parse_data_final.py : Prétraitement des spectres des molécules d'un fichier mgf, pour pouvoir effectuer différentes analyses.
* visualisation_clusters.py : Permet une visualisation (plot ou web) des molécules dans les clusters (sans prendre en compte les distances inter et intra cluster).
* visualisation_mol_communes.py : Permet de visualiser l'intersection de 2 clustering en regardant les clusters ayant des molécules en communs.
* visualisation_similarite.py : Permet la visualisation des similarités en utilisant 3 algorithmes: UMAP, T_SNE et PCA.


## Modules nécessaires

Plusieurs dépendances sont nécessaires au fonctionnement du programme. Elles sont contenus dans le fichier `requirements.txt`.
Pour les installer, il est possible d'utiliser la commande suivante:
```bash
% pip install -r requirements.txt
```

## Utilisation du programme

Pour lancer le programme vous pouvez utiliser la commande suivante:
```bash
% python main.py <command>
```

Par exemple
```bash
% python main.py cosinus
```

L'attribut `<command>` peut prendre les valeurs suivantes:
* `parse` : Permet de faire le parsing des fichiers `Cluster.mgf` et `ALL_GNPS_cleaned.mgf`.
* `cosinus` : Calcule les matrices de similarités cosinus sur les spectres.
* `fingerprint` : Calcule les matrices de similarités Tanimoto sur les fingerprints des SMILES.
* `groupes` : Calcule les matrices de similarités Tanimoto sur les fingerprints des SMILES.
* `mcl` : Fait le clustering de Markov Clustering sur les 3 types de similarités.
* `hdbscan` : Fait le clustering HDBSCAN sur les 3 types de similarités.

L'exécution doit être fait dans un certain ordre.
Pour pouvoir calculer les similarités (commandes `cosinus`, `fingerprint` et `groupes`) il faut avoir les fichiers contenant les molécules et donc avoir exécuté le parsing avec `% python main.py parse`.
Pour utiliser les méthodes de clustering (commandes `mcl` et `hdbscan`) il faut avoir des fichiers de similarités.
Pour effectuer la visualisation (commandes `cluster` et `intersection`), il faut des fichiers de clustering.


Pour la visualisation, il y a 2 commandes possibles qui prennent chacunes des arguments:
* `cluster` : Permet de visualiser les molécules contenues dans les clusters. Il a besoin d'un fichier issue du clustering au format .txt.
```bash
python3 main.py cluster -p <cluster>
```
* `intersection` : Permet de visualiser l'intersection des clusters à partir de deux fichiers de clustering. Il a besoin de deux fichiers issues du clustering au format .txt.
```bash
python3 main.py cluster -p <cluster1> -s <cluster2>
```

Exemple d'utilisation:
```bash
python3 main.py cluster -p cluster_molecules/HDBSCAN_cosinus_spectre/energy_37.0_precursor_M+Na.txt

python3 main.py intersection -p cluster_molecules/HDBSCAN_cosinus_spectre/energy_37.0_precursor_M+Na.txt -s cluster_molecules/HDBSCAN_fingerprints_smiles/energy_37.0_precursor_M+Na.txt
```


Les différentes commandes sont exécutés sur des fichiers mgf que nous avons choisis. La liste est la suivante:
```bash
energy_50.0_precursor_M+H.mgf
energy_10.0_precursor_M+H.mgf
energy_30.0_precursor_M+H.mgf
energy_37.0_precursor_M+Na.mgf
```
