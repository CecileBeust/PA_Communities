# import modules
import pandas as pd
import os
import argparse
import sys
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import scipy.stats as stats
from functions_enrichment import extract_seeds, load_geneage, extract_genes_from_cluster, extract_genes_from_comm, remove_seeds, create_mapping_file, getMappingDict, map_ensembl_to_symbol, fisher, hypergeome

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

from utilities import create_dico_disease_seeds, build_communities_list, filter_cluster, create_cluster_dico, load_networks

# Argparse
parser = argparse.ArgumentParser(
    prog="enrichment_physio_aging.py", 
    description="functions to perform enrichment of gene communities with physiological aging genes from the GenAge database and lists of genes up and down-regulated during physiological aging"
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
print(len(list_ids_analyzed))

all_nodes = load_networks(comm_path)
seeds = list(extract_seeds(orpha_codes, list_ids_analyzed))
genage = load_geneage('Data_PhysioAging/genage_human.csv', seeds, all_nodes)

dico_cluster_diseases = create_cluster_dico(cluster_output)
filtered_dico_cluster = filter_cluster(dico_cluster_diseases)

dico_comm_nodes = extract_genes_from_comm(comm_path, 100, list_ids_analyzed)
dico_clusters_nodes = extract_genes_from_cluster(comm_path, filtered_dico_cluster, 100)
print(dico_clusters_nodes.keys())

# mapping of gene identifiers
mapping_file_path = create_mapping_file(path + '/Data_PhysioAging/GSE103232_hs_blood_batch2_counts_rpkm.xls')

def create_enrichment_file_genage(genage: list, dico_clusters_nodes: dict, background: int) -> None:
    """Function who perform Fisher and Hypergeometric tests to assess the enrichment of 
    physiological aging genes from the GenAge database in the 6 clusters of gene communities

    Args:
        genage (list): preprocessed list of physiological aging genes from GenAge
        dico_clusters_nodes (dict): dictionray of the diseases-associated communities in each cluster
        background (int): number of unique nodes in the multiplex network, taken as backgroun for
        the statistical significance of the analysis
    """
    df = pd.DataFrame(np.zeros((7, 3)))
    df.columns = ['Cluster', 'Fisher test p-value', 'Hypergeometric test p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        print(cluster)
        nodes_cluster = dico_clusters_nodes[cluster]
        set_gene_age = set(genage)
        overlap = len(set_gene_age.intersection(set(nodes_cluster)))
        cluster_unique = len(nodes_cluster) - overlap
        deg_unique = len(set_gene_age) - overlap                                       
        print("HYPERGEOMETRIC TEST")
        h_pval = hypergeome(nodes_cluster, genage, background)
        print("FISCHER TEST")
        f_pval = fisher(nodes_cluster, genage, background)
        df._set_value(i, 'Cluster', int(cluster[8:]))
        print(df)
        df._set_value(i, 'Fisher test p-value', f_pval)
        df._set_value(i, 'Hypergeometric test p-value', h_pval)
        i += 1
    df.to_csv(path + f"output_tables/Enrichment_GeneAge.tsv", sep="\t", index=False)
    print(df)

create_enrichment_file_genage(genage=genage, dico_clusters_nodes=dico_clusters_nodes, background=len(all_nodes))

def extract_list_genes_DEG_physio_aging(file: str) -> tuple[list, list, list]:
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
    print(file)
    print("Up-regulated genes:", len(up_genes))
    print("Down-regulated genes:", len(down_genes))
    print("Other genes :", len(other_genes))
    return up_genes, down_genes, other_genes

def create_filtered_enrichment_lists_physio_aging_DEG(file: str, seeds_list: list, is_up: int, all_nodes: list) -> list:
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
    genes_up, genes_down, other_genes = extract_list_genes_DEG_physio_aging(file)

    # If we want to analyze UP regulated genes
    if is_up == 1:
        # map to gene symbol
        genes_up_GS = map_ensembl_to_symbol(mapping_file_path, genes_up, 'ID', 'external_gene_id')
        print(f"{len(genes_up_GS)} genes UP")
        # remove seeds
        gene_up_wo_seeds = remove_seeds(genes_up_GS, seeds_list)
        print(f"{len(gene_up_wo_seeds)} genes UP without the seeds")
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
        genes_down_GS = map_ensembl_to_symbol(mapping_file_path, genes_down, 'ID', 'external_gene_id')
        print(f"{len(genes_down_GS)} genes DOWN")
        # remove seeds
        gene_down_wo_seeds = remove_seeds(genes_down_GS, seeds_list)
        print(f"{len(gene_down_wo_seeds)} genes DOWN without seeds")
        # check if aging genes are present in networks
        genes_down_wo_seeds_in_ntw = []
        for gene in gene_down_wo_seeds:
            if gene in all_nodes:
                genes_down_wo_seeds_in_ntw.append(gene)
        print(f"{len(genes_down_wo_seeds_in_ntw)} genes DOWN present in multiplex (WO seeds)")
        return genes_down_wo_seeds_in_ntw



def create_enrichment_files(deg: str, tissue: str, background: int):
    if deg == "UP":
        genes_enrich = create_filtered_enrichment_lists_physio_aging_DEG(
            file=f'human-{tissue}.txt',
            seeds_list=seeds,
            is_up=1,
            all_nodes=all_nodes
            )
    elif deg == "DOWN":
        genes_enrich = create_filtered_enrichment_lists_physio_aging_DEG(
            file=f'human-{tissue}.txt',
            seeds_list=seeds,
            is_up=0,
            all_nodes=all_nodes
            )

    df = pd.DataFrame(np.zeros((7, 3)))
    df.columns = ['Cluster', 'Fisher test p-value', 'Hypergeometric test p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        print(" ")
        print(cluster)
        nodes_cluster = set(dico_clusters_nodes[cluster])
        set_DEG = set(genes_enrich)
        overlap = len(set_DEG.intersection(nodes_cluster))
        cluster_unique = len(nodes_cluster) - overlap
        deg_unique = len(set_DEG) - overlap
        
        print("####### HYPERGEOMETRIC TEST #########")
        h_pval = hypergeome(nodes_cluster, set_DEG, background)
        print("######## FISHER #############")
        f_pval = fisher(nodes_cluster, set_DEG, background)
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'Fisher test p-value'] = f_pval
        df.at[i, 'Hypergeometric test p-value'] = h_pval
        i += 1
    df.to_csv(path + f"output_tables/Enrichment_genes_{deg}_{tissue}.tsv", sep="\t", index=False)
    print(df)

def enrich_clusters_all_physio_aging(background: int, cluster: str, matrix: np.array, genage: list, tissue="GenAge", deg="GenAge"):
    if deg == "up":
        genes_enrich = create_filtered_enrichment_lists_physio_aging_DEG(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                seeds_list=seeds,
                is_up=1,
                all_nodes=all_nodes
                )
    elif deg == "down":
        genes_enrich = create_filtered_enrichment_lists_physio_aging_DEG(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                seeds_list=seeds,
                is_up=0,
                all_nodes=all_nodes
                )
    elif deg == "GenAge" or tissue == "GenAge":
        genes_enrich = genage
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
    
    if tissue == "GenAge" and deg == "GenAge":
        j = 0
    elif tissue == "blood" and deg == "up":
        j = 1
    elif tissue == "blood" and deg == "down":
        j = 2
    elif tissue == "skin" and deg == "up":
        j = 3
    elif tissue == "skin" and deg == "down":
        j = 4
    elif tissue == "brain" and deg == "up":
        j = 5
    elif tissue == "brain" and deg == "down":
        j = 6
    elif tissue == "muscle" and deg == "up":
        j = 7
    elif tissue == "muscle" and deg == "down":
        j = 8
    elif tissue == "breast" and deg == "up":
        j = 9
    elif tissue == "breast" and deg == "down":
        j = 10
    h_pval = hypergeome(dico_clusters_nodes[cluster], genes_enrich, background)
    matrix[i][j] = hypergeome(dico_clusters_nodes[cluster], genes_enrich, background)
    return h_pval


def heatmap_enrichment(dico_cluster_nodes, background):
    enrichment_matrix = np.zeros((len(dico_clusters_nodes), 11))
    df = pd.DataFrame(columns = ['Cluster', 'GenAge', 'Blood up', 'Blood down', 'Skin up', 'Skin down', 'Brain up', 'Brain down', 'Muscle up', 'Muscle down', 'Breast up', 'Breast down'])
    i = 0
    for cluster in dico_clusters_nodes:
        pval_geneage = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "GenAge", "GenAge")
        pval_blood_up = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "blood", "up")
        pval_blood_down = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "blood", "down")
        pval_skin_up = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "skin", "up")
        pval_skin_down = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "skin", "down")
        pval_brain_up = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "brain", "up")
        pval_brain_down = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "brain", "down")
        pval_muscle_up = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "muscle", "up")
        pval_muscle_down = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "muscle", "down")
        pval_breast_up = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "breast", "up")
        pval_breast_down = enrich_clusters_all_physio_aging(background, cluster, enrichment_matrix, genage, "breast", "down")
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'GenAge'] = pval_geneage
        df.at[i, 'Blood up'] = pval_blood_up
        df.at[i, 'Blood down'] = pval_blood_down
        df.at[i, 'Skin up'] = pval_skin_up
        df.at[i, 'Skin down'] = pval_skin_down
        df.at[i, 'Brain up'] = pval_brain_up
        df.at[i, 'Brain down'] = pval_brain_down
        df.at[i, 'Muscle up'] = pval_muscle_up
        df.at[i, 'Muscle down'] = pval_muscle_down
        df.at[i, 'Breast up'] = pval_breast_up
        df.at[i, 'Breast down'] = pval_breast_down
        i += 1
    df.to_csv(f"output_tables/Enrichment_clusters_physio_aging_genes.tsv", sep="\t", index=False)
    print(enrichment_matrix)
    ax = plt.axes()
    y_labels = ["Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"]
    x_labels = [
        "GenAge",
        "Blood \n up-reg. genes", 
        "Blood \n down-reg. genes", 
        "Skin \n up-reg. genes",
        "Skin \n  down-reg. genes",
        "Brain \n up-reg. genes",
        "Brain \n down-reg. genes",
        "Muscle \n up-reg. genes",
        "Muscle \n down-reg. genes",
        "Breast \n up-reg. genes",
        "Breast \n down-reg. genes"
        ]
    plt.tick_params(axis='both', which='major', labelsize=8, labelbottom = False, bottom=False, top = False, labeltop=True)
    vmin = 0
    vmax = 0.05
    mask = np.logical_or(enrichment_matrix < vmin, enrichment_matrix > vmax)
    heatmap = sns.heatmap(
        enrichment_matrix,
        cmap="crest_r",
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