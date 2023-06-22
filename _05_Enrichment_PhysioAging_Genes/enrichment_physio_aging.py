# import modules
import pandas as pd
import os
import argparse
import sys
import networkx as nx
import numpy as np
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import seaborn as sns
import scipy.stats as stats
from functions_enrichment import load_networks, extract_seeds, load_geneage, extract_genes_from_cluster, extract_genes_from_comm, fisher, hypergeome

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, build_communities_list, filter_cluster, create_cluster_dico

# Argparse
parser = argparse.ArgumentParser(
    prog="enrichment_physio_aging.py", 
    description="functions to perform enrichment of gene communities with physiological aging genes from the GenAge database"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# define path to data
data_folder = os.path.join(os.path.dirname(__file__), '..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')

(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
#list_id_analyzed = ['93932', '920', '909', '904', '90348', '90324', '90322', '90321', '902', '90154', '90153', '897', '895', '894', '808', '79477', '79476', '79474', '79333', '79325', '79087', '769', '758', '75496', '740', '633', '500', '498359', '487825', '477', '435628', '412057', '363618', '357074', '3455', '3437', '33445', '33364', '3322', '3163', '3051', '2963', '287', '286', '2834', '280679', '280365', '2658', '263534', '263487', '2500', '220295', '2078', '2067', '1942', '1901', '1860', '1807', '1775', '163746', '1387', '137608', '1340', '1299', '1297', '101028', '100']
print(len(list_ids_analyzed))

all_nodes = load_networks(comm_path)
seeds = list(extract_seeds(orpha_codes, list_ids_analyzed))
geneage = load_geneage('Data_PhysioAging/genage_human.csv', seeds, all_nodes)

dico_cluster_diseases = create_cluster_dico(cluster_output)
filtered_dico_cluster = filter_cluster(dico_cluster_diseases)

dico_comm_nodes = extract_genes_from_comm(comm_path, 100, list_ids_analyzed)
dico_clusters_nodes = extract_genes_from_cluster(comm_path, filtered_dico_cluster, 100)
print(dico_clusters_nodes.keys())

######################################
####### GENEAGE ######################
######################################

def make_enrichment_geneage(geneage: list, dico_clusters_nodes: dict, background: int) -> None:
    """Function who perform Fisher and Hypergeometric tests to assess the enrichment of 
    physiological aging genes from the GenAge database in the 6 clusters of gene communities

    Args:
        geneage (list): preprocessed list of physiological aging genes from GenAge
        dico_clusters_nodes (dict): dictionray of the diseases-associated communities in each cluster
        background (int): number of unique nodes in the multiplex network, taken as backgroun for
        the statistical significance of the analysis
    """
    df = pd.DataFrame(np.zeros((7, 3)))
    df.columns = ['Cluster', 'Fisher p-value', 'Hypergeometric p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        print(" ")
        print(cluster)
        nodes_cluster = dico_clusters_nodes[cluster]
        set_gene_age = set(geneage)
        overlap = len(set_gene_age.intersection(set(nodes_cluster)))
        cluster_unique = len(nodes_cluster) - overlap
        deg_unique = len(set_gene_age) - overlap
        venn2(
            subsets=(cluster_unique, deg_unique, overlap), 
            set_labels=(f"{cluster}", "GeneAge genes"),
            set_colors=('red', 'green')
            )                                              
        plt.show()
        print("HYPERGEOMETRIC TEST")
        h_pval = hypergeome(nodes_cluster, geneage, background)
        print("FISCHER TEST")
        f_pval = fisher(nodes_cluster, geneage, background)
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'Fisher p-value'] = f_pval
        df.at[i, 'Hypergeometric p-value'] = h_pval
        i += 1
    df.to_csv(path + f"output_tables/Enrichment_GeneAge.tsv", sep="\t", index=False)
    print(df)

make_enrichment_geneage(geneage=geneage, dico_clusters_nodes=dico_clusters_nodes, background=len(all_nodes))

def create_mapping_file(rpkm_file: str):
    """Function to create a mapping file (Ensembl gene symbols to HUGO
    gene symbols) from the RPKM file of the study of Irizar et. al

    Args:
        rpkm_file (str): path to the RPKM file
    """
    df = pd.read_excel(rpkm_file)
    df = df.iloc[:, :2]
    print(df)
    df.to_csv(path + "/Data_PhysioAging/Mapping_Ensembl_GeneSymbol.txt", sep="\t", index=None)

create_mapping_file(path + '/Data_PhysioAging/GSE103232_hs_blood_batch2_counts_rpkm.xls')
mapping_file_path = path + '/Data_PhysioAging/Mapping_Ensembl_GeneSymbol.txt'

def getMappingDict(filePath, convertFrom, convertTo):
	df=pd.read_csv(filePath, sep="\t")
	df=df[[convertFrom, convertTo]]
	df=df.dropna(axis=0)
	if('ID' in df.columns):
		df['ID']= df['ID'].astype(str)
		df['ID']= df['ID'].astype(str)
	
	df=df.set_index(convertFrom)
	
	mappingDict=df.to_dict()[convertTo]
	
	return mappingDict

def extract_list_genes(file: str) -> tuple[list, list, list]:
    """Function to extract the lists of genes up, down or up and down-regulated
    during physiological aging described in the study of Irizar et. al

    Args:
        enrichment_file (str): path to the file containing the list of genes

    Returns:
        tuple[list, list, list]: list of genes up-regulated during physiological 
        aging in the studied tissue, list of genes down-regulated, and list of 
        genes both up and down-regulated. 
    """
    up_genes = []
    down_genes = []
    other_genes = []
    is_up = None
    with open(file, 'r') as f:
        lines = f.read().splitlines()

        for line in lines:
            line = line.strip()

            if line == "Up:":
                is_up = True
            elif line == "Down:":
                is_up = False
            elif line == "Other:":
                 is_up = None
            elif is_up is True:
                up_genes.append(line)
            elif is_up is False:
                down_genes.append(line)
            else:
                other_genes.append(line)
    print("Up-regulated genes:", len(up_genes))
    print("Down-regulated genes:", len(down_genes))
    print("Other genes :", len(other_genes))
    return up_genes, down_genes, other_genes

def map_ensembl_to_symbol(list_genes: list, From: str, To: str) -> list:
    """Function to map Ensembl genes identifiers to HUGO gene symbols

    Args:
        list_genes (list): genes to map
        From (str): identifiers to be mapped (Ensembl)
        To (str): identifiers to map (HUGO gene symbols)

    Returns:
        list: list of genes identified with HUGO genes symbols
    """
    mapping_dict = getMappingDict(
    filePath=mapping_file_path,
    convertFrom=From,
    convertTo=To
    )
    new_list = []
    for gene in list_genes:
        if gene in mapping_dict:
            new_list.append(mapping_dict[gene])
    return new_list

def remove_seeds(gene_list: list, seeds_list: list) -> list:
    """Function which remove seeds nodes used to identify the communities
    from a list of genes

    Args:
        gene_list (list): list of genes to analyze
        seeds_list (list): list of seeds to remove

    Returns:
        list: new list of genes without seeds
    """
    new_list = []
    for gene in gene_list:
        if not gene in seeds_list:
            new_list.append(gene)
    return new_list

def create_enrichment_lists(file: str, seeds_list: list, is_up: int, all_nodes: list) -> list:
    """Function to create the lists of genes up, down or both up and down regulated
    during physiological aging for the enrichments analysis. The seeds are removed from
    the list. 

    Args:
        file (str): path to file containing the list of genes
        seeds_list (list): list of seeds nodes
        is_up (int): 1 if we perform enrichment with genes up-regulated, 0 if we 
        performe enrichment with genes down-regulated
        all_nodes (list): list of all nodes in the multiplex network

    Returns:
        list: list of genes for the enrichment analysis
    """
    # extract lists of genes from file
    genes_up, genes_down, other_genes = extract_list_genes(file)

    # If we want to analyze UP regulated genes
    if is_up == 1:
        # map to gene symbol
        genes_up_GS = map_ensembl_to_symbol(genes_up, 'ID', 'external_gene_id')
        print(len(genes_up_GS))
        # remove seeds
        gene_up_wo_seeds = remove_seeds(genes_up_GS, seeds_list)
        print(len(gene_up_wo_seeds))
        # check if aging genes are present in networks
        genes_up_wo_seeds_in_ntw = []
        for gene in gene_up_wo_seeds:
            if gene in all_nodes:
                genes_up_wo_seeds_in_ntw.append(gene)
        print(f"{len(genes_up_wo_seeds_in_ntw)} genes UP present in multiplex (WO seeds)")
        return genes_up_wo_seeds_in_ntw
    
    # If we want to analyze DOWN regulated genes
    elif is_up == 0:
        # map to gene names
        genes_down_GS = map_ensembl_to_symbol(genes_down, 'ID', 'external_gene_id')
        print(len(genes_down_GS))
        # remove seeds
        gene_down_wo_seeds = remove_seeds(genes_down_GS, seeds_list)
        print(len(gene_down_wo_seeds))
        # check if aging genes are present in networks
        genes_down_wo_seeds_in_ntw = []
        for gene in gene_down_wo_seeds:
            if gene in all_nodes:
                genes_down_wo_seeds_in_ntw.append(gene)
        print(f"{len(genes_down_wo_seeds_in_ntw)} genes DOWN present in multiplex (WO seeds)")
        return genes_down_wo_seeds_in_ntw
    
    # If we want to analyze both UP and DOWN regulated genes
    else:
        # add genes up and down
        whole_degs = list(set(genes_up + genes_down))
        print(len(whole_degs))
        # map to gene symbols
        whole_degs_GS = map_ensembl_to_symbol(whole_degs, 'ID', 'external_gene_id')
        print(len(whole_degs_GS))
        # remove seeds
        whole_degs_wo_seeds = remove_seeds(whole_degs_GS, seeds_list)
        # check if aging genes are present in networks
        whole_degs_wo_seeds_in_ntw = []
        for gene in whole_degs_wo_seeds:
            if gene in all_nodes:
                whole_degs_wo_seeds_in_ntw.append(gene)
        print(f"{len(whole_degs_wo_seeds_in_ntw)} genes UP or DOWN present in multiplex (WO seeds)")
        return whole_degs_wo_seeds_in_ntw


def create_enrichment_files(deg: str, tissue: str, background: int):
    if deg == "UP":
        genes_enrich = create_enrichment_lists(
            file=f'human-{tissue}.txt',
            seeds_list=seeds,
            is_up=1,
            all_nodes=all_nodes
            )
        venn_color = 'green'
    elif deg == "DOWN":
        genes_enrich = create_enrichment_lists(
            file=f'human-{tissue}.txt',
            seeds_list=seeds,
            is_up=0,
            all_nodes=all_nodes
            )
        venn_color = 'tomato'
    elif deg == "UP+DOWN":
        genes_enrich = create_enrichment_lists(
            file=f'human-{tissue}.txt',
            seeds_list=seeds,
            is_up=2,
            all_nodes=all_nodes
            )
        venn_color = 'dodgerblue'

    df = pd.DataFrame(np.zeros((7, 3)))
    df.columns = ['Cluster', 'Fisher p-value', 'Hypergeometric p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        print(" ")
        print(cluster)
        nodes_cluster = set(dico_clusters_nodes[cluster])
        set_DEG = set(genes_enrich)
        overlap = len(set_DEG.intersection(nodes_cluster))
        cluster_unique = len(nodes_cluster) - overlap
        deg_unique = len(set_DEG) - overlap
        venn2(
            subsets=(cluster_unique, deg_unique, overlap), 
            set_labels=(f"Cluster {cluster}", f"Genes {deg} in {tissue}"),
            set_colors=('red', venn_color)
            )                  
        plt.show()
        
        print("####### HYPERGEOMETRIC TEST #########")
        h_pval = hypergeome(nodes_cluster, set_DEG, background)
        print("######## FISHER #############")
        f_pval = fisher(nodes_cluster, set_DEG, background)
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'Fisher p-value'] = f_pval
        df.at[i, 'Hypergeometric p-value'] = h_pval
        i += 1
    df.to_csv(path + f"output_tables/Enrichment_genes_{deg}_{tissue}.tsv", sep="\t", index=False)
    print(df)

#create_enrichment_files("DOWN", "skin", all_nodes)
#create_enrichment_files("DOWN", "blood", all_nodes)
#create_enrichment_files("WHOLE", "blood", all_nodes)

def enrich_cluster(background: int, cluster, matrix, tissue="GeneAge", deg="GeneAge"):
    if deg == "up":
        genes_enrich = create_enrichment_lists(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                seeds_list=seeds,
                is_up=1,
                all_nodes=all_nodes
                )
    elif deg == "down":
        genes_enrich = create_enrichment_lists(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                seeds_list=seeds,
                is_up=0,
                all_nodes=all_nodes
                )
    elif deg == "up+down":
        genes_enrich = create_enrichment_lists(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                seeds_list=seeds,
                is_up=2,
                all_nodes=all_nodes
                )
    elif deg == "GeneAge" or tissue == "GeneAge":
        genes_enrich = geneage
    if cluster == "cluster_1":
        i = 0
    elif cluster == "cluster_2":
        i = 1
    elif cluster == "cluster_3":
        i = 2
    elif cluster == "cluster_4":
        i = 3
    elif cluster == "cluster_5":
        i = 4
    elif cluster == "cluster_6":
        i = 5
    
    if tissue == "GeneAge" and deg == "GeneAge":
        j = 0
    elif tissue == "blood" and deg == "up":
        j = 1
    elif tissue == "blood" and deg == "down":
        j = 2
    elif tissue == "blood" and deg == "up+down":
        j = 3
    elif tissue == "skin" and deg == "up":
        j = 4
    elif tissue == "skin" and deg == "down":
        j = 5
    elif tissue == "skin" and deg == "up+down":
        j = 6
    h_pval = hypergeome(dico_clusters_nodes[cluster], genes_enrich, background)
    matrix[i][j] = hypergeome(dico_clusters_nodes[cluster], genes_enrich, background)
    return h_pval


def heatmap_enrichment(dico_cluster_nodes, background):
    enrichment_matrix = np.zeros((len(dico_clusters_nodes), 7))
    df = pd.DataFrame(columns = ['Cluster', 'GeneAge', 'Blood up', 'Blood down', 'Blood whole', 'Skin up', 'Skin down', 'Skin whole'])
    i = 0
    for cluster in dico_clusters_nodes:
        pval_geneage = enrich_cluster(background, cluster, enrichment_matrix, "GeneAge", "GeneAge")
        pval_blood_up = enrich_cluster(background, cluster, enrichment_matrix, "blood", "up")
        pval_blood_down = enrich_cluster(background, cluster, enrichment_matrix, "blood", "down")
        pval_blood_whole = enrich_cluster(background, cluster, enrichment_matrix, "blood", "up+down")
        pval_skin_up = enrich_cluster(background, cluster, enrichment_matrix, "skin", "up")
        pval_skin_down = enrich_cluster(background, cluster, enrichment_matrix, "skin", "down")
        pval_skin_whole = enrich_cluster(background, cluster, enrichment_matrix, "skin", "up+down")
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'GeneAge'] = pval_geneage
        df.at[i, 'Blood up'] = pval_blood_up
        df.at[i, 'Blood down'] = pval_blood_down
        df.at[i, 'Blood whole'] = pval_blood_whole
        df.at[i, 'Skin up'] = pval_skin_up
        df.at[i, 'Skin down'] = pval_skin_down
        df.at[i, 'Skin whole'] = pval_skin_whole
        i += 1
    df.to_csv(f"Enrichment_clusters_aging_genes.tsv", sep="\t", index=False)
    print(enrichment_matrix)
    ax = plt.axes()
    y_labels = ["Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"]
    x_labels = [
        "GeneAge",
        "Genes UP \n BLOOD", 
        "Genes DOWN \n BLOOD", 
        "Genes \n UP+DOWN  \n BLOOD",
        "Genes UP \n SKIN",
        "Gene DOWN \n SKIN",
        "Genes \n UP+DOWN \n SKIN"
        ]
    plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)
    vmin = 0
    vmax = 0.05
    mask = np.logical_or(enrichment_matrix < vmin, enrichment_matrix > vmax)
    heatmap = sns.heatmap(
        enrichment_matrix,
        cmap=sns.cubehelix_palette(light = .8, n_colors=10, as_cmap = False, reverse = True),
        mask=mask,
        vmin=vmin,
        vmax=vmax,
        xticklabels=x_labels,
        yticklabels=y_labels,
        ax=ax
    )
    plt.savefig(path + "output_figures/Heatmap_PhysioAging.png", bbox_inches='tight')
    plt.show()

heatmap_enrichment(dico_cluster_nodes=dico_clusters_nodes, background=len(all_nodes))