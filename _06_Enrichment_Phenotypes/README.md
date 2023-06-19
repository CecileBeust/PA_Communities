# Enrichment of clusters using HPO phenotypes

This folder contains Python codes to enrich clusters of diseases using phenotypes from the Human Phenotype Ontology.

## Files

* ```HP.csv``` : Ontology file of HPO terms
* ```PA_diseases.tsv``` : File containing the Premature Aging diseases analyzed
* ```clusters_100.tsv``` : File describing the diseases in each cluster
* ```enrichment_phenotypes.py``` : Python script containing enrichment functions. The enrichment results are filtered using orsum. 

## Usage

Download the phenotypes file corresponding the diseases analyzed from HPO and store it in a folder named ```HPOTermsOrphaDiseases ```. Then run the script ```enrichment_phenotypes.py```. 

## Output

* ```EnrichPhenoResults```: Folder containing the enrichment results for HPO phenotypes in each cluster of diseases.
* Orsum folders containing files listing signicant phenotypes in each clusters
* ```ME_results_phenotypes```: Folder containing the results of the filtereing of enrichment results with orsum. Contains the heatmaps of most significant terms.

## References

* Köhler, S. et al. The Human Phenotype Ontology in 2017. en. Nucleic Acids Research 45,
D865–D876. ISSN: 0305-1048, 1362-4962. doi:10.1093/nar/gkw1039 (Jan. 2017)

* Ozisik, O., Térézol, M. & Baudot, A. orsum: a Python package for filtering and comparing
enrichment analyses using a simple principle. en. BMC Bioinformatics 23, 293. ISSN:
1471-2105. doi:10.1186/s12859-022-04828-2 (Dec. 2022).