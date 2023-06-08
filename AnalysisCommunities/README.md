# Analysis, Clustering, Enrichment and Visualization of gene communities

This folder contains different scripts to analyze communities. The scripts can be used individually depending on the type of analysis you want to perform on you communities, or can be used altogther in the same pipeline.

## Folders

* ```EnrichmentAgingGenes```: This folder contains codes to perform an enrichment of the PA associated communities in physiological aging genes. The code can be used and adapted to assess the enrichment of other genes sets in communities. 

* ```EnrichmentPhenotypes```: This folder containg codes to perform an enrichment of diseases in clusters using HPO phentoypes

* ```PlotCommunities```: This folder contains codes to generate files for the visualisation of communities in the Cytoscape application.

* ```data```: Data folder containing the list of PA diseases and their associated genes from ORPHANET
 
* ```gmt_files```: Folder containing GMT files for the enrichment analysis with g:Profiler 

## Files

* ```analysis_final_communities.py```: Python functions that are used in other scripts
* ```analysis_genes_in_clusters.py```: contain functions to analyze the genes present in each cluster of communities
* ```analysis_genes_in_communities.py```: contain functions to analyze the genes present in the communities
* ```cluster_communities.py```: contain functions to cluster communities of genes based on Jaccard index
* ```enrichment_clusters.py```: contain functions to enrich the clusters using biological annotations
* ```enrichment_communities.py```: contain functions to enrich the communities using biogical annotations
* ```orsum_enrichment_clusters.py```: contain functions to filter the enrichment results in biological annotations for the clusters


## References

* Köhler, S. et al. The Human Phenotype Ontology in 2017. en. Nucleic Acids Research 45,
D865–D876. ISSN: 0305-1048, 1362-4962. doi:10.1093/nar/gkw1039 (Jan. 2017).

* Rath, A. et al. Representation of rare diseases in health information systems: The orphanet approach to serve a wide range of end users. Human Mutation 33, 803–808. ISSN:10597794. doi:10.1002/humu.22078 (2012).

* Raudvere, U. et al. g:Profiler: a web server for functional enrichment analysis and con-
versions of gene lists (2019 update). en. Nucleic Acids Research 47, W191–W198. ISSN:
0305-1048, 1362-4962. doi:10.1093/nar/gkz369 (July 2019)

* Shannon, P. et al. Cytoscape : A Software Environment for Integrated Models of Biomolec-
ular Interaction Networks. Genome Research 1, 2498–2504. doi:10.1101/gr.1239303.
metabolite (2003).