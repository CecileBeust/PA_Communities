# Enrichment of gene communities, or clusters of communities using biological annotations 

This folder contains codes to perform an enrichment of gene communities (or clusters of communities) using biological annotations from the Gene Ontology (GO Biological Processes and GO Cellular Components) and Reactome pathways terms. It is possible to do the same analysis using other biological annotations.

The enrichement analysis are performed using g:Profiler, and the enrichement results are filtered using orsum.

## Files 

* ```enrichment_clusters.py```: Python script containing functions to perform enrichment analysis of clusters of communities. 
* ```enrichment_communities.py```: Python script containing functions to perform enrichment analysis of individual communities
* ```orsum_enrichment_clusters.py```: Python script containing functions to filter the enrichment results using orsum.

## Usage

```python enrichment_clusters.py -p /path/where/communities/folders/are/stored```

```python enrichment_communities.py -p /path/where/communities/folders/are/stored```

After having modified the script with your version of orsum (line 53), run:

```python orsum_enrichment_clusters.py```


## Output

The script ```enrichment_clusters.py``` generates enrichment results file for each cluster of communities (files ```cluster_100_ID.tsv``` where ID is the identifier of the cluster). It uses the file ```cluster_output_100_0.7.tsv``` generated after the clustering of communities, which is available in the ```_00_data``` folder. 

The script ```enrichment_communities.py``` generates enrichment results file for each community, and a file ```enrichment_communities_100.xlsx``` that recapitulates the enrichment of all the communities.

The script ```orsum_enrichment_clusters.py``` will create files for the analysis with orsum, for each cluster. Also it will generate folders containing the results of the multiple enrichment analysis, with the corresponding heatmaps.

## References

* Ashburner, M. et al. Gene Ontology: tool for the unification of biology. en. Nature Genet-
ics 25. Number: 1 Publisher: Nature Publishing Group, 25–29. ISSN: 1546-1718. doi:10.
1038/75556 (May 2000).

* The Gene Ontology Consortium et al. The Gene Ontology knowledgebase in 2023. Ge-
netics 224, iyad031. ISSN: 1943-2631. doi:10.1093/genetics/iyad031 (May 2023).

* Raudvere, U. et al. g:Profiler: a web server for functional enrichment analysis and con-
versions of gene lists (2019 update). en. Nucleic Acids Research 47, W191–W198. ISSN:
0305-1048, 1362-4962. doi:10.1093/nar/gkz369 (July 2019).

* Ozisik, O., Térézol, M. & Baudot, A. orsum: a Python package for filtering and comparing
enrichment analyses using a simple principle. en. BMC Bioinformatics 23, 293. ISSN:
1471-2105. doi:10.1186/s12859-022-04828-2 (Dec. 2022).