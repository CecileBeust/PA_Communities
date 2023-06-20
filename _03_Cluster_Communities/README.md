# Clustering of gene communities based on Jaccard index

This folder contains the Python codes allowing to perform a clustering of gene communities based on their similarity using Jaccard index.

The Jaccard index is the ratio of the common genes between two communities over the union of their genes.

## Files

* ```cluster_communities.py```: Python script containing functions to compute pairwise Jaccard indices between communities and perform the clustering of communities.

* ```analyse_genes_in_clusters.py```: Python script containing functions to analyse the genes represented inside clusters of communities.

## Usage

```python cluster_communities.py -p /path/where/communities/folders/are/stored```

```python analyse_genes_in_clusters.py -p /path/where/communities/folders/are/stored```

## Output

The script ```cluster_communities.py``` will generate a clustering of communities based on Jaccard index. It is possible to download to corresponding dendrogram and clustermap.

The script ```analyse_genes_in_clusters.py``` will generate an excel file of the genes represented in each cluster of communities.