# Identification, clustering and enrichment of communities in a multiplex network: application to Premature Aging diseases 

This repository presents a pipeline allowing to:
* Identify communities associated to Premature Aging (PA) diseases in a multiplex biological network  
* Perform a clustering of the communities based on Jaccard index
* Perform enrichment analysis of the clusters using biological annotation or other lists of genes or phenoypes

Here we identifiy communities related to Premature Aging diseases, using their associated genes from ORPHANET as seeds in a multiplex network to perform an iterative Random Walk with Restart (itRWR). 

The itRWR algorithm is available as a Python package on GitHub : https://github.com/anthbapt/itRWR

This pipeline can be used and/or adapted to any study requiring the identification of communities in a multiplex network and their analysis.

## Folders

* ```_00_data```: Folder containing the data files needed for the analysis (diseases identifiers files, gmt files and some output files used in different steps of the analysis)
* ```_01_Community_Identification```: Folder containing codes for the community identification with itRWR. It contains a subset toy example of a multiplex biological network composed of 4 layers : Protein-protein interactions, Molecular complexes, Pathways and Coexpression networks. This subset is extracted from complete networks, which are available on the NDEx server: https://www.ndexbio.org/index.html#/search?searchType=All&searchString=cecile.beust&searchTermExpansion=false

* ```_02_Analysis_Communities```: Folder containing codes to analyze communities and visualize them in Cytoscape

* ```_03_Cluster_Communities```: Folder containing the codes to perform a clustering of communities based on Jaccard index.

* ```_04_Enrichment_BioAnnotations```: Folder containing the codes to perform enrichent analysis of the communities, or clusters of communities, using biological annotations.

* ```_05_Enrichment_PhysioAging_Genes```: Folder containing the codes to perfom enrichment analysis of the clusters using physiological aging genes.

* ```_06_Enrichment_Phenotypes```: Folder containing codes to perform enrich;ent analysus of the clusters using phenotypes

* ```Tests``` : Folder containing a test file to check the results of the community identification algorithm.

## Usage

The codes provided here are adapted to the study of Premature Aging diseases identified with their ORPHANET identifiers. It is possible to adapt the codes for other application cases. 

## References

* Rath, A. et al. Representation of rare diseases in health information systems: The or-
phanet approach to serve a wide range of end users. Human Mutation 33, 803–808. ISSN:10597794. doi:10.1002/humu.22078 (2012).