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
import useful_functions_enrichment

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

# variables statement
(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]

all_nodes = list(load_networks(comm_path))
seeds = list(useful_functions_enrichment.extract_seeds(orpha_codes, list_ids_analyzed))
genage = useful_functions_enrichment.load_geneage('Data_PhysioAging/genage_human.csv', seeds, all_nodes)

dico_cluster_diseases = create_cluster_dico(cluster_output)
filtered_dico_cluster = filter_cluster(dico_cluster_diseases)

dico_comm_nodes = useful_functions_enrichment.extract_genes_from_comm(comm_path, 100, list_ids_analyzed)
dico_clusters_nodes = useful_functions_enrichment.extract_genes_from_cluster(comm_path, filtered_dico_cluster, 100)

# mapping of gene identifiers in list of DEG physio aging genes
mapping_file_path = useful_functions_enrichment.create_mapping_file(path + '/Data_PhysioAging/GSE103232_hs_blood_batch2_counts_rpkm.xls')

def create_enrichment_file_genage(genage: list, dico_clusters_nodes: dict, all_nodes: list) -> None:
    """Function who performs Fisher and Hypergeometric tests to assess the enrichment of 
    physiological aging genes from the GenAge database in the 6 clusters of gene communities

    Args:
        genage (list): preprocessed list of physiological aging genes from GenAge
        dico_clusters_nodes (dict): dictionray of the diseases-associated communities in each cluster
        background (int): number of unique nodes in the multiplex network, taken as background for
        the statistical significance of the analysis
    """
    df = pd.DataFrame(np.zeros((6, 3)))
    df.columns = ['Cluster', 'Fisher test p-value', 'Hypergeometric test p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        nodes_cluster = dico_clusters_nodes[cluster]                         
        h_pval = useful_functions_enrichment.hypergeome(list1=nodes_cluster, list2=genage, gene_pool=all_nodes)
        f_pval = useful_functions_enrichment.fisher(list1=nodes_cluster, list2=genage, gene_pool=all_nodes)
        df._set_value(i, 'Cluster', int(cluster[8:]))
        df._set_value(i, 'Fisher test p-value', f_pval)
        df._set_value(i, 'Hypergeometric test p-value', h_pval)
        i += 1
    df.to_csv(path + f"output_tables/enrichment_clusters_GenAge.csv", sep=",", index=False)

create_enrichment_file_genage(genage=genage, dico_clusters_nodes=dico_clusters_nodes, all_nodes=all_nodes)

def enrich_clusters_all_physio_aging(dico_clusters_nodes: dict, all_nodes: list, seeds: list, cluster: str, matrix: np.array, genage: list, mapping_file_path: str, tissue="GenAge", deg="GenAge"):
    """Function that performs an hypergeometric test for the enrichment analysis of a list of
    physiological aging genes in the 6 clusters

    Args:
        dico_clusters_nodes (dict): dico containing the nodes of each cluster
        all_nodes (list): list of unique nodes in the multiplex network
        seeds (list): list of seeds nodes
        cluster (str): identifier of the cluster
        matrix (np.array): enrichment matrix
        genage (list): list of physiological aging genes from GenAge
        mapping_file_path (str): path to the mapping file used to map Ensembl gene identifiers to HUGO gene symbols
        tissue (str, optional): name of the tissue for the enrichment of DEG genes. Defaults to "GenAge".
        deg (str, optional): "up" or "down". Defaults to "GenAge".

    Returns:
        float: p-value of the hypergeometric test
    """
    if deg == "up":
        genes_enrich = useful_functions_enrichment.create_filtered_enrichment_lists_physio_aging_DEG(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                mapping_file_path=mapping_file_path,
                seeds_list=seeds,
                is_up=1,
                all_nodes=all_nodes
                )
    elif deg == "down":
        genes_enrich = useful_functions_enrichment.create_filtered_enrichment_lists_physio_aging_DEG(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                mapping_file_path=mapping_file_path,
                seeds_list=seeds,
                is_up=0,
                all_nodes=all_nodes
                )
    elif deg == "GenAge" or tissue == "GenAge":
        genes_enrich = genage
    
    # state coordinates of enrichment results in the matrix
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
    # extract nodes of cluster
    nodes_cluster = dico_clusters_nodes[cluster]
    # perform hypergeometric test
    h_pval = useful_functions_enrichment.hypergeome(list1=nodes_cluster, list2=genes_enrich, gene_pool=all_nodes)
    # add p-value in the matrix
    matrix[i][j] = h_pval
    return h_pval


def heatmap_enrichment(dico_clusters_nodes: dict, all_nodes: list, seeds: list, genage: list, mapping_file_path: str):
    # initialize enrichment matrix
    enrichment_matrix = np.zeros((len(dico_clusters_nodes), 11))
    df = pd.DataFrame(columns = ['Cluster', 'GenAge', 'Blood up-reg. genes', 'Blood down-reg. genes', 'Skin up-reg. genes', 'Skin down-reg. genes', 'Brain up-reg. genes', 'Brain down-reg. genes', 'Muscle up-reg. genes', 'Muscle down-reg. genes', 'Breast up-reg. genes', 'Breast down-reg. genes'])
    i = 0
    for cluster in dico_clusters_nodes:
        # perform hypergeometric tests for all lists of genes
        pval_genage = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "GenAge", "GenAge")
        pval_blood_up = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "blood", "up")
        pval_blood_down = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "blood", "down")
        pval_skin_up = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "skin", "up")
        pval_skin_down = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "skin", "down")
        pval_brain_up = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "brain", "up")
        pval_brain_down = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "brain", "down")
        pval_muscle_up = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "muscle", "up")
        pval_muscle_down = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "muscle", "down")
        pval_breast_up = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "breast", "up")
        pval_breast_down = enrich_clusters_all_physio_aging(dico_clusters_nodes, all_nodes, seeds, cluster, enrichment_matrix, genage, mapping_file_path, "breast", "down")
        df.at[i, 'Cluster'] = str(cluster)
        df.at[i, 'GenAge'] = pval_genage
        df.at[i, 'Blood up-reg. genes'] = pval_blood_up
        df.at[i, 'Blood down-reg. genes'] = pval_blood_down
        df.at[i, 'Skin up-reg. genes'] = pval_skin_up
        df.at[i, 'Skin down-reg. genes'] = pval_skin_down
        df.at[i, 'Brain up-reg. genes'] = pval_brain_up
        df.at[i, 'Brain down-reg. genes'] = pval_brain_down
        df.at[i, 'Muscle up-reg. genes'] = pval_muscle_up
        df.at[i, 'Muscle down-reg. genes'] = pval_muscle_down
        df.at[i, 'Breast up-reg. genes'] = pval_breast_up
        df.at[i, 'Breast down-reg. genes'] = pval_breast_down
        i += 1
    # export enrichment results to table
    df.to_csv(f"output_tables/enrichment_clusters_physio_aging_genes.csv", sep=",", index=False)
    # plot enrichment heatmap
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

heatmap_enrichment(dico_clusters_nodes=dico_clusters_nodes, all_nodes=all_nodes, seeds=seeds, genage=genage, mapping_file_path=mapping_file_path)