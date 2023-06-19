from analyse_final_communities import create_dico_disease_seeds
from cluster_communities import build_communities_list, build_similarity_matrix
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
import argparse
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# Argparse
parser = argparse.ArgumentParser(
    prog="analyse_genes_in_communities.py", 
    description="functions to analyse the gene composition of communities"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")

(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "orpha_codes_PA.txt")
(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
similarity_matrix_10 = build_similarity_matrix(communities_10)[0]

def analyse_genes_in_comm(communities: list, list_ids_analyzed: list) -> tuple[list, list, dict]:
    """
    Function to get the list of total genes in a set of communities, the list of genes without the seeds and a 
    dictionary of genes and their belonging communities

    Args:
        communities (list) : list of paths to the communities to analyze
        list_ids_analyzed (list) : list of the ORPHANET identifiers of diseases communities to analyze

    Return:
        (tuple) : the list of total genes in the communities, the list of genes without seeds, the dicionary of the genes
        and their belonging communities
    """
    dico_gene_comm = {}
    genes_total = []
    genes_wo_seeds = []
    genes_in_several_comm = []
    seeds = set()
    df = pd.read_table(path + "orpha_codes_PA.txt")
    for index, row in df.iterrows():
        if row[0] in list_ids_analyzed:
            seeds_disease = row[1].split(",")
            for seed in seeds_disease:
                seeds.add(seed)
    for comm, id in zip(communities, list_ids_analyzed):
        with open(comm, 'r') as file:
            for line in file:
                gene = line.rstrip()
                # if the gene has not been registered in the lits of total genes
                if not gene in genes_total:
                    # add it to the list of genes and register it in the dico
                    genes_total.append(gene)
                    dico_gene_comm[gene] = [comm]
                # if the gene has been registered
                else:
                    dico_gene_comm[gene] += [comm]
                    # if the gene has not been registered in the list of genes appearing in several modules
                    if gene not in genes_in_several_comm:
                        # add it to the list of genes involved in several modules
                        genes_in_several_comm.append(gene)
                # if the gene is not a seed
                if not gene in seeds:
                    if not gene in genes_wo_seeds:
                        genes_wo_seeds.append(gene)
    return genes_total, genes_wo_seeds, dico_gene_comm

(genes_total, genes_wo_seeds, dico_gene_comm) = analyse_genes_in_comm(communities_10, [740, 902])
print(f"Total genes in commmunities : {genes_total}")
print(f"Genes without seeds in communities : {genes_wo_seeds}")
print(f"Dico genes in communities : {dico_gene_comm}")

def generate_table_genes_in_comm(dico_gene_comm, genes_wo_seeds) -> None:
    """Function to generate a table of the genes in a set of 
    communities and the number of communities they belong to

    Args:
        dico_gene_comm (dict): dictionary of genes and their
        associated communities
        genes_wo_seeds (list): list of genes in the communities
        without the seeds
    Return:
        None
    """
    df_100 = pd.DataFrame(columns=['Gene name', 'Number of PA communities memberships'])
    i = 0
    for gene in dico_gene_comm.keys():
        if gene in genes_wo_seeds:
            df_100._set_value(i, 'Gene name', gene)
            nb_comm = len(dico_gene_comm[gene])
            df_100._set_value(i, 'Number of PA communities memberships', nb_comm)
            i += 1
    df_100.sort_values(by=['Number of PA communities memberships'], ascending=False)
    print(df_100)
    #df_100.to_csv(path + "genes_comm.tsv", sep="\t", index=None)

generate_table_genes_in_comm(dico_gene_comm, genes_wo_seeds)

def generate_excel_genes(dico_gene_comm: dict) -> None:
    """Function to generate an excel file gathering 
    information about the genes inside a set of communities
    (number of communities it belongs to, degrees in the different networks of 
    the multiplex, maximum degree, sum of the degree). The function also prints 
    the results of a Pearson correlation test to assess the correlation between
    the number of communities a gene belongs to and its degrees in the
    networks

    Args:
        dico_gene_comm (dict): dictionary of genes and their associated communities
    """
    df = pd.DataFrame(columns=['Gene symbol', 'Nb of communities', 'degree PPI ntw', 'degree Pathways ntw', 'degree Coexp ntw', 'degree Complexes ntw', 'max degree', 'sum degrees'])
    i = 0
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/Pathways.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes.tsv", create_using = nx.Graph)
    for gene in dico_gene_comm:
        df._set_value(i, 'Gene symbol', gene)
        nb_comm = 0
        comm_names = []
        for comm in dico_gene_comm[gene]:
            nb_comm += 1
            comm_names.append(comm)
            df._set_value(i, 'Nb of communities', nb_comm)
        all_degs = []
        for elt in ppi.degree([gene]):
            deg_ppi = elt[1]
            all_degs.append(deg_ppi)
            df._set_value(i, 'degree PPI ntw', deg_ppi)
        for elt in pathways.degree([gene]):
            deg_pathways = elt[1]
            all_degs.append(deg_pathways)
            df._set_value(i, 'degree Pathways ntw', deg_pathways)
        for elt in coexp.degree([gene]):
            deg_coexp = elt[1]
            all_degs.append(deg_coexp)
            df._set_value(i, 'degree Coexp ntw', deg_coexp)
        for elt in complexes.degree([gene]):
            deg_complexes = elt[1]
            all_degs.append(deg_complexes)
            df._set_value(i, 'degree Complexes ntw', deg_complexes)
        df._set_value(i, 'max degree', max(all_degs))
        df._set_value(i, 'sum degrees', sum(all_degs))
        i += 1
    print(df)
    # Pearson correlation test
    comm = df['Nb of communities']
    max_deg = df['max degree']
    sum_deg = df['sum degrees']
    corr_max, p_value_max = pearsonr(comm, max_deg)
    print('Correlation Nb comm / max deg:', corr_max)
    print('P-value:', p_value_max)
    corr_sum, p_value_sum = pearsonr(comm, sum_deg)
    print('Correlation Nb comm / sum deg:', corr_sum)
    print('P-value:', p_value_sum)
    df.to_csv("gene_in_communities.tsv", sep="\t")

generate_excel_genes(dico_gene_comm)