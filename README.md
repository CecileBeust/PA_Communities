# Identification, clustering and enrichment of communities in a multiplex network : application to Premature Aging diseases 

This repository contains codes to :
* Identify communities associated to Premature Aging (PA) diseases in a multiplex biological network  
* Perform a clustering of the communities based on Jaccard index
* Perform enrichment analysis of the clusters using biological annotation or other lists of genes or phenoypes

Here we identifiy communities related to Premature Aging diseases, using their associated genes as seeds in a multiplex network to perform an iterative random walk with restart (itRWR). The itRWR algorithm is available as a Python package on GitHub : https://github.com/anthbapt/itRWR

This pipeline can be used and/or adapted to any study requiring the identification of communities in a multiplex network and their analysis.

## Folders

* ```AnalysisCommunities``` : Folder containing codes to analyze communities. Codes for the clustering of the communities based on Jaccard Index, and the enrichment analysis of the obtained clusters are provided, as well as codes for the visualisation of communities in Cytoscape.
* ```IDCommunity``` : Folder containing codes for the community identification with itRWR. It contains a toy example of a multiplex network composed of 4 layers : Protein-protein interactions, Molecular complexes, Pathways and Coexpression networks. 
* ```Tests``` : Folder containing a test file to check the results of the community identification algorithm.
