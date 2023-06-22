# Enrichment of clusters of gene communities using physiological aging genes

This folder contains codes to enrich clusters of gene communities using lists of physiological aging genes. Here the lists of physiological aging genes used are the following : 

* Genes from the GenAge database
* Genes up and down-regulated during physiological aging in blood and skin, from the study of Irizar et al.

## Folders

* ```Data_PhysioAging```: Folder containing the physiological aging data.

This folder contains the following data files : 

* ```GSE103232_hs_blood_batch2_counts_rpkm.xls```: File containing Ensembl identifiers of genes to map to HUGO genes symbols, from the study of Irizar et al. 
* ```geneage_human.csv```: physiological aging genes extracted from GenAge
* ```human-blood.txt```: list of genes differentially expressed during physiological aging in blood (from the study of Irizar et al)
* ```human-skin.txt```: list of genes differentially expressed during physiological aging in skin (from the study of Irizar et al)
* ```Mapping_Ensembl_GeneSymbol.txt```: mapping file for Ensembl gene symbols tp HUGO gene symbols

## Files

* ```enrichment_physio_aging.py```: Python script containing the functions to perform the enrichment analysis
* ```functions_enrichment.py```: Python script containing other functions used in ```enrichment_physio_aging.py```

## Usage

```python enrichment_physio_aging.py -p /path/where/communities/folders/are/stored```

## Output

* ```Enrichment_clusters_aging_genes.tsv``` : file containing the p-values of hypergeometric tests for the enrichment of the cluster in physiological aging genes

The code also allows to create a heatmap of enrichment results. 

## References

* Aramillo Irizar, P. et al. Transcriptomic alterations during ageing reflect the shift from cancer to degenerative diseases in the elderly. en. Nature Communications 9. Number: 1 Publisher: Nature Publishing Group, 327. ISSN: 2041-1723. doi:10.1038/s41467-017-
02395-2 (Jan. 2018).

* Tacutu, R. et al. Human Ageing Genomic Resources: new and updated databases. en.
Nucleic Acids Research 46, D1083–D1090. ISSN: 0305-1048, 1362-4962. doi:10.1093/
nar/gkx1042 (Jan. 2018)