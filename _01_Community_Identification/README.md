# Identification of communities associated with Premature Aging diseases in a multiplex network

This folder contains the Python scripts allowing running the identification of network communities associated with the 67 Premature Aging (PA) diseases extracted from ORPHANET. 

The communities are identified using an iterative Random Walk with Restart algorithm, itRWR (https://github.com/anthbapt/itRWR.git), which is based on the MultiXrank Python package. The itRWR algorithm is applied to a multiplex network of biological interactions using the PA disease causative genes as seeds. 

## Data

* Folder ```multiplex```: contains a multiplex network of biological interactions composed of 4 layers (PPI, Pathways, Co-expression, Complexes).

* ```orpha_codes_PA.txt```: file containing the Premature Aging diseases and their associated causative genes. Data were extracted from HPO and ORPHANET. We removed from the file the causative genes that are not present in the biological multiplex network. Indeed, these causative genes cannot be used as seeds in the itRWR. 

## Files

* ```config_ID.yml```: configuration files for each PA disease analyzed, with ID being ORPHANET codes of the diseases.  

* ```seeds_ID.txt```: files containing the causative genes associated with each disease. These causative genes are used as seeds by the itRWR algorithm.


* ```run_CI.py```: Python script allowing running the identification of communities for the 67 PA diseases associated with at least one causative gene (i.e., associated with at least one seed) in the multiplex network. Here we use a number of iterations of 100, but changing this parameter in the script is possible. This will change the size of the obtained community.

## Usage

    python run_CI.py

## Output

For each disease, a folder named ```results_100_ID``` will be created, with ID being the ORPHANET code of the disease. This folder contains the following files:

* ```config.yml```: copy of the configuration file used for the disease
* ```multiplex_1.tsv```: a file containing the rankings of all the nodes of the multiplex network after running itRWR using the causative gene(s) associated with the disease as seed(s)
* ```seeds_ID.txt```: a file containing the nodes of the community identified for the disease. 

## References

* Baptista, A., Gonzalez, A. & Baudot, A. Universal multilayer network exploration by ran-
dom walk with restart. en. Communications Physics 5, 170. ISSN: 2399-3650. doi:10 .
1038/s42005-022-00937-9 (July 2022).

* Rath, A. et al. Representation of rare diseases in health information systems: The or-
phanet approach to serve a wide range of end users. Human Mutation 33, 803–808. ISSN: 10597794. doi:10.1002/humu.22078 (2012).

* Köhler S, et al. The Human Phenotype Ontology in 2021. Nucleic Acids Res. 2021 Jan 8;49(D1):D1207-D1217. doi: 10.1093/nar/gkaa1043. PMID: 33264411; PMCID: PMC7778952.
