# Identification of communities associated to Premature Aging diseases in a multiplex network

This folder contains the Python script allowing to run the identification of communities associated to 67 Premature Aging diseases from ORPHANET in a multiplex network of biological interactions. 

The communities are identified using and iterative Random Walk with Restart algorithm, itRWR (https://github.com/anthbapt/itRWR.git), based on the MultiXrank Python package.

## Files

* ```config_ID.yml```: configuration files for each disease analyzed, with ID being ORPHANET codes of the diseases.  

* ```seeds_ID.txt```: files containing the causative genes associated to each diseases. These causative genes are used as seeds by the itRWR algorithm.

* ```orpha_codes_PA.txt```: file containing the Premature Aging diseases analyzed and their associated causative genes from ORPHANET. We removed from the file the causative genes which were nor present in our multiplex network of biological interactions, and thus can not be used as seeds. 

* ```run_CI.py```: Python script allowing to run the identification of communities for the 67 ORPHANET diseases associated to at least one causative gene (one seed) in the multiplex network. Here we use a number of iteration of 100, but it is possible to change this parameter in the script. This will change the size of the obtained communities.

## Usage

```python run_CI.py```

## Output

For each disease, a folder named ```results_100_ID``` will be created, with ID being the ORPHANET code of the disease. This folder contains the following files:

* ```config.yml```: copy of the configuration file used for the disease
* ```multiplex_1.tsv```: a file containing the rankings of each node of the multiplex network after itRWR was run using the causative gene(s) of the disease as seed(s)
* ```seeds_ID.txt```: a file containing the nodes of the community identified for the disease. 

## References

* Baptista, A., Gonzalez, A. & Baudot, A. Universal multilayer network exploration by ran-
dom walk with restart. en. Communications Physics 5, 170. ISSN: 2399-3650. doi:10 .
1038/s42005-022-00937-9 (July 2022).

* Rath, A. et al. Representation of rare diseases in health information systems: The or-
phanet approach to serve a wide range of end users. Human Mutation 33, 803–808. ISSN: 10597794. doi:10.1002/humu.22078 (2012).