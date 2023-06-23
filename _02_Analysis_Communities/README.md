# Analysis and visualization of communities 

This folder contains different scripts to analyze communities of genes, here identified for Premature Aging diseases in a multiplex network of biological interaction using an iterative Random Walk with Restart algorithm (itRWR)


## Files

* ```analysis_genes_in_communities.py```: contains Python functions to analyze the genes present in the communities
* ```visualization_communities.py```: contains Python functions to generate tabulated files for the visualization of communities in Cytoscape

## Usage

    python analysis_genes_in_communities.py -p /path/where/communities/folders/are/stored

    python visualization_communities.py -p /path/where/communities/folders/are/stored

## Output

The script ```analysis_genes_in_communities``` generates a table of the genes inside communities in the ```output_tables``` folder. 
The script ```visualization_communities``` generates tsv files of communities, allowing them to be visualized individually or as group of communities in Cytoscape. 

<img src="visualization_3comm.png" alt="Alt text" title="Communities of Hutchinson-Gilford Progeria Syndrome (ORHPANET code 740), Ataxia telangiectasia (ORPHANET code 100) and Classical Ehlers-Danlos Syndrome (ORPHANET code 287) visualized in Cytoscape">


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