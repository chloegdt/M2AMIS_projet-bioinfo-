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
* parse_data_final.py : Prépare et organise les spectres moléculaires pour les différentes analyses.
* visualisation_clusters.py : [Description du script]
* visualisation_mol_communes.py : [Description du script]
* visualisation_similarite.py : [Description du script]


## Modules nécessaires

Plusieurs dépendances sont nécessaires au fonctionnement du programme. Elles sont contenus dans le fichier `requirements.txt`.
Pour les installer, il est possible d'utiliser la commande:
```bash
% pip install -r requirements.txt
```

## Utilisation du programme

Pour lancer le programme vous pouvez utiliser la commande suivante:
```bash
% python main.py <command>
```

L'attribut `<command>` peut prendre les valeurs suivantes:
* `parse` : Permet de faire le parsing des fichiers `Cluster.mgf` et `ALL_GNPS_cleaned.mgf`.
* `cosinus` : Calcule les matrices de similarités cosinus sur les spectres.
* `fingerprint` : Calcule les matrices de similarités Tanimoto sur les fingerprints des SMILES.
* `groupes` : Calcule les matrices de similarités Tanimoto sur les fingerprints des SMILES.
* `mcl` : Fait le clustering de Markov Clustering sur les 3 types de similarités.
* `hdbscan` : Fait le clustering HDBSCAN sur les 3 types de similarités.
