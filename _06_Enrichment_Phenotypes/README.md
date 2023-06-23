# Enrichment of clusters using HPO phenotypes

This folder contains Python codes to enrich clusters of diseases using phenotypes from the Human Phenotype Ontology.

## Data

* The ```Data_HPO``` folder contains the HPO ontology file ```HP.csv``` and the phenotypes files from HPO associated to the 67 Premature Aging diseases analyzed.

## Files

* ```enrichment_phenotypes.py``` : Python script containing enrichment functions. The enrichment results are filtered using orsum. 
* ```hpo_phenotypes```: GMT file of HPO phenotypes and associated diseases, built with the ```enrichment_phenotypes.py``` script. 

## Usage

```python enrichment_phenotypes.py ```

## Output

* ```output_tables```: Folder containing the enrichment results for HPO phenotypes in each cluster of diseases.
* ```outout_orsum```: Folder containing the results of the filtereing of enrichment results with orsum. Contains the heatmaps of most significant terms.

## References

* Köhler, S. et al. The Human Phenotype Ontology in 2017. en. Nucleic Acids Research 45,
D865–D876. ISSN: 0305-1048, 1362-4962. doi:10.1093/nar/gkw1039 (Jan. 2017)

* Ozisik, O., Térézol, M. & Baudot, A. orsum: a Python package for filtering and comparing
enrichment analyses using a simple principle. en. BMC Bioinformatics 23, 293. ISSN:
1471-2105. doi:10.1186/s12859-022-04828-2 (Dec. 2022).