# Identification, clustering and enrichment of communities in a multiplex network : application to Premature Aging diseases 

This repository contains codes to identify, cluster and enrich premature aging disease-associated communities in a multiplex biological network. 

Here we identifiy communities related to Premature Aging (PA) diseases, using their associated genes as seeds in a multiplex network to perform an iterative random walk with restart (itRWR).

This pipeline can be used for any study related to the identification of communities of genes in a multiplex network, their analysis, clustering or functional enrichment.

## Folders

* ```AnalysisCommunities``` : Folder containing codes to analyze communities. Codes for the clustering of the communities based on Jaccard Index, and the enrichment analysis of the obtained clusters are provided, as well as codes for the visualisation of communities in Cytoscape.
* ```IDCommunity``` : Folder containing codes for the community identification with itRWR (cf https://github.com/anthbapt/itRWR). It contains a toy example of a multiplex network composed of 4 layers : Protein-protein interactions, Molecular complexes, Pathways and Coexpression networks. 
*```Tests``` : Folder containing a test file to check the results of the community identification algorithm.
