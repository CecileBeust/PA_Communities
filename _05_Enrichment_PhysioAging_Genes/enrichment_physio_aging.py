# import modules
import pandas as pd
import os
import argparse
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
from statsmodels.sandbox.stats.multicomp import multipletests
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
all_seeds = list(useful_functions_enrichment.extract_seeds(orpha_codes, list_ids_analyzed))
# MODIF KEEP SEEDS
# genage = useful_functions_enrichment.load_geneage('Data_PhysioAging/genage_human.csv', seeds, all_nodes)
genage = useful_functions_enrichment.load_geneage_keep_seeds('Data_PhysioAging/genage_human.csv', all_seeds, all_nodes)

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
    df.columns = ['Cluster', 'Hypergeometric test p-value', 'Corrected p-value']
    i = 0
    for cluster in dico_clusters_nodes:
        nodes_cluster = dico_clusters_nodes[cluster]                         
        h_pval = useful_functions_enrichment.hypergeome(list1=nodes_cluster, list2=genage, gene_pool=all_nodes)
        df._set_value(i, 'Cluster', int(cluster[8:]))
        df._set_value(i, 'Hypergeometric test p-value', h_pval)
        i += 1
    pvals = df['Hypergeometric test p-value'].to_list()
    p_adjusted = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
    df['Corrected p-value'] = p_adjusted
    df.to_csv(path + f"output_tables/enrichment_clusters_GenAge.csv", sep=",", index=False)

create_enrichment_file_genage(genage=genage, dico_clusters_nodes=dico_clusters_nodes, all_nodes=all_nodes)

def enrich_clusters_all_physio_aging(dico_clusters_nodes: dict, all_nodes: list, seeds: list, cluster: str, genage: list, mapping_file_path: str, tissue="GenAge", deg="GenAge"):
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
        genes_enrich = useful_functions_enrichment.create_filtered_enrichment_lists_physio_aging_DEG_keep_seeds(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                mapping_file_path=mapping_file_path,
                seeds_list=seeds,
                is_up=1,
                all_nodes=all_nodes
                )
    elif deg == "down":
        genes_enrich = useful_functions_enrichment.create_filtered_enrichment_lists_physio_aging_DEG_keep_seeds(
                file=f'Data_PhysioAging/human-{tissue}.txt',
                mapping_file_path=mapping_file_path,
                seeds_list=seeds,
                is_up=0,
                all_nodes=all_nodes
                )
    elif deg == "GenAge" or tissue == "GenAge":
        genes_enrich = genage

    # extract nodes of cluster
    nodes_cluster = dico_clusters_nodes[cluster]
    intersection = set(nodes_cluster).intersection(set(genes_enrich))

    # perform hypergeometric test
    h_pval = useful_functions_enrichment.hypergeome(list1=nodes_cluster, list2=genes_enrich, gene_pool=all_nodes)

    return h_pval, nodes_cluster, genes_enrich, intersection

def fill_enrichment_table(df: pd.DataFrame, position: int, cluster: str, physio_aging_genes: str, pval: float, nodes_cluster: list, genes_enrich: list, intersection: set, seeds: list):
    df.at[position, 'Cluster'] = cluster
    df.at[position, "List of physio aging genes"] = physio_aging_genes
    df.at[position, 'Nb genes cluster'] = len(list(nodes_cluster))
    df.at[position, 'Nb genes in physio aging list'] = len(list(genes_enrich))
    df.at[position, "Background"] = len(all_nodes)
    df.at[position, "Nb genes intersection"] = len(intersection)
    if intersection != set():
        df.at[position, "Intersection"] = intersection
    df.at[position, "p-value"] = pval
    physio_aging_seeds = list()
    for gene in genes_enrich:
        if gene in seeds:
            physio_aging_seeds.append(gene)
    df.at[position, "Seeds"] = physio_aging_seeds

def heatmap_enrichment(dico_clusters_nodes: dict, all_nodes: list, seeds: list, genage: list, mapping_file_path: str):
    physio_aging_genes_list = ['GenAge', 'Blood up-reg. genes', 'Blood down-reg. genes', 'Skin up-reg. genes', 'Skin down-reg. genes', 'Brain up-reg. genes', 'Brain down-reg. genes', 'Muscle up-reg. genes', 'Muscle down-reg. genes', 'Breast up-reg. genes', 'Breast down-reg. genes']
    tissues_list = ['GenAge', 'blood', 'blood', 'skin', 'skin', 'brain', 'brain', 'muscle', 'muscle', 'breast', 'breast']
    deg_list = ['GenAge', 'up', 'down', 'up', 'down', 'up', 'down', 'up', 'down', 'up', 'down']
    summary_pvalues = pd.DataFrame(columns = ['Cluster', 'GenAge', 'Blood up-reg. genes', 'Blood down-reg. genes', 'Skin up-reg. genes', 'Skin down-reg. genes', 'Brain up-reg. genes', 'Brain down-reg. genes', 'Muscle up-reg. genes', 'Muscle down-reg. genes', 'Breast up-reg. genes', 'Breast down-reg. genes'])
    i = 0
    for cluster in dico_clusters_nodes:
        dico_cluster_diseases = create_cluster_dico(cluster_output)
        filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
        seeds_cluster = useful_functions_enrichment.select_seeds_from_cluster(dico_disease_seeds=dico_disease_seeds, dico_cluster=filtered_dico_cluster, cluster_id=cluster)
        enrich_results_table = pd.DataFrame(columns=["Cluster", "List of physio aging genes", "Nb genes cluster", "Nb genes in physio aging list", "p-value", "Corrected p-value", "Background", "Nb genes intersection", "Intersection", "Seeds"])
        summary_pvalues.at[i, 'Cluster'] = str(cluster)
        j = 1
        for physio_aging_genes, tissue, deg in zip(physio_aging_genes_list, tissues_list, deg_list):
            print(physio_aging_genes, tissue, deg)
            # perform hypergeometric tests for all lists of genes
            pval, nodes_cluster, genes_enrich, intersection = enrich_clusters_all_physio_aging(dico_clusters_nodes=dico_clusters_nodes, all_nodes=all_nodes, seeds=seeds_cluster, cluster=cluster, genage=genage, mapping_file_path=mapping_file_path, tissue=tissue, deg=deg)
            fill_enrichment_table(df=enrich_results_table, position=j, cluster=cluster, physio_aging_genes=physio_aging_genes, pval=pval, nodes_cluster=nodes_cluster, genes_enrich=genes_enrich, intersection=intersection, seeds=seeds_cluster)
            j += 1
        # adjust the pvalues: Benjamini-Hochberg
        pvals = enrich_results_table['p-value'].to_list()
        p_adjusted = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        enrich_results_table['Corrected p-value'] = p_adjusted
        # export enrichment file for the cluster
        enrich_results_table.to_csv(f"output_tables/enrichment_details/enrichment_{cluster}.csv", sep=",", index=False)

        # fill p-values summary table for the cluster
        summary_pvalues.at[i, 'GenAge'] = p_adjusted[0]
        summary_pvalues.at[i, 'Blood up-reg. genes'] = p_adjusted[1]
        summary_pvalues.at[i, 'Blood down-reg. genes'] = p_adjusted[2]
        summary_pvalues.at[i, 'Skin up-reg. genes'] = p_adjusted[3]
        summary_pvalues.at[i, 'Skin down-reg. genes'] = p_adjusted[4]
        summary_pvalues.at[i, 'Brain up-reg. genes'] = p_adjusted[5]
        summary_pvalues.at[i, 'Brain down-reg. genes'] = p_adjusted[6]
        summary_pvalues.at[i, 'Muscle up-reg. genes'] = p_adjusted[7]
        summary_pvalues.at[i, 'Muscle down-reg. genes'] = p_adjusted[8]
        summary_pvalues.at[i, 'Breast up-reg. genes'] = p_adjusted[9]
        summary_pvalues.at[i, 'Breast down-reg. genes'] = p_adjusted[10]
        i += 1
    # export p-values summary table to csv file
    summary_pvalues.to_csv(f"output_tables/enrichment_clusters_physio_aging_genes.csv", sep=",", index=False)

    # create an enrichment matrix for heatmap visualization
    summary_pvalues = summary_pvalues.drop(columns=["Cluster"])
    #enrichment_matrix = np.zeros((len(dico_clusters_nodes), 11))
    enrichment_matrix = np.zeros(summary_pvalues.shape)
    enrichment_matrix[enrichment_matrix == 0] = summary_pvalues.values[enrichment_matrix == 0]

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
    plt.tick_params(axis='both', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True)
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

heatmap_enrichment(dico_clusters_nodes=dico_clusters_nodes, all_nodes=all_nodes, seeds=all_seeds, genage=genage, mapping_file_path=mapping_file_path)

def clean_and_merge_enrichment_physio_aging_files(dico_clusters_nodes: dict):
    csv_dir = Path(path + "output_tables/enrichment_details/")
    csv_data = {}
    for csv_file in csv_dir.glob('*.csv'):
        csv_name = csv_file.stem
        csv_data[csv_name] = pd.read_csv(csv_file, sep=",")
    writer = pd.ExcelWriter(path + f"output_tables/enrichment_details/enrichment_physioaging.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in csv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.close()
    for cluster in dico_clusters_nodes:
        os.remove(f"output_tables/enrichment_details/enrichment_{cluster}.csv")

clean_and_merge_enrichment_physio_aging_files(dico_clusters_nodes=dico_clusters_nodes)

def create_number_DEG_genes_supp_file():
    list_files = ['Data_PhysioAging/human-blood.txt', 'Data_PhysioAging/human-skin.txt', 'Data_PhysioAging/human-brain.txt', 'Data_PhysioAging/human-muscle.txt', 'Data_PhysioAging/human-breast.txt']
    list_tissues = ['Blood', 'Skin', 'Brain', 'Muscle', 'Breast']
    supp_file = pd.DataFrame(columns=["Tissue", "Number of up-regulated genes", "Number of down-regulated genes", "Number of 'other' genes"])
    i = 0
    for file, tissue in zip(list_files, list_tissues):
        genes_up, genes_down, other_genes = useful_functions_enrichment.extract_list_genes_DEG_physio_aging(file)
        supp_file._set_value(i, "Tissue", tissue)
        supp_file._set_value(i, "Number of up-regulated genes", len(genes_up))
        supp_file._set_value(i, "Number of down-regulated genes", len(genes_down))
        supp_file._set_value(i, "Number of 'other' genes", len(other_genes))
        i += 1
    print(supp_file)
    supp_file.to_csv(path + f"output_tables/number_of_DEG_genes_Irizar_et_al.csv", sep=",", header=1, index=False)

create_number_DEG_genes_supp_file()