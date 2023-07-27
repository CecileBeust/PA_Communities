"""
Functions to analyze the gene composition of communities
and generate associated files
"""

# import modules
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
import argparse
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
import sys
from pathlib import Path

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, build_communities_list

# define path to data
data_folder = os.path.join(os.path.dirname(__file__), '..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')


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

# variables statement
(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)
print(f"Dico diseases seeds: {dico_disease_seeds}")
print(" ")
print(f"List all diseases: {list_id}")
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
print(" ")
print(f"Number of communities: {len(communities_100)}")
print(" ")
print(f"Diseases not analyzed: {not_analyzed}")
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]

def analyse_genes_in_comm(communities: list, list_ids_analyzed: list) -> tuple[list, list, dict]:
    """
    Function to get the list of total genes in a set of communities, the list of genes without the seeds and a 
    dictionary of genes and their belonging communities

    Args:
        communities (list) : list of paths to the communities to analyze
        list_ids_analyzed (list) : list of the ORPHANET identifiers of diseases communities analyzed

    Return:
        (tuple) : the list of total genes in the communities, the list of genes without seeds, the dicionary of the genes
        and their belonging communities
    """
    dico_gene_comm = {}
    genes_total = []
    genes_wo_seeds = []
    genes_in_several_comm = []
    seeds = set()
    df = pd.read_csv(orpha_codes, sep="\t", header=None)
    for index, row in df.iterrows():
        if str(row[0]) in list_ids_analyzed:
            seeds_disease = row[1].split(",")
            for seed in seeds_disease:
                seeds.add(seed)
    # check if there are as much communities as diseases analyzed
    assert len(communities) == len(list_ids_analyzed)
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

(genes_total, genes_wo_seeds, dico_gene_comm) = analyse_genes_in_comm(communities_100, list_ids_analyzed)
print(" ")
print(f"Total genes in commmunities : {len(genes_total)}")
print(" ")
print(f"Genes without seeds in communities : {len(genes_wo_seeds)}")

def generate_excel_genes(dico_gene_comm: dict) -> None:
    """Function to generate file gathering 
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
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
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
    df_sorted = df.sort_values(by=['Nb of communities'], ascending=False)
    df_sorted.to_csv("output_tables/genes_in_communities.csv", sep=",")
    
    
generate_excel_genes(dico_gene_comm)

def list_genes_in_communities(orpha_file: str, list_comm: list, list_ids_analyzed: list, comm_path: str):
    diseases_names = pd.read_csv(orpha_file, sep="\t", header=None)
    diseases = list()
    for index, row in diseases_names.iterrows():
        if row[0][6:] in list_ids_analyzed:
            diseases.append(row[1])
    for comm, disease, id in zip(list_comm, diseases, list_ids_analyzed):
        df = pd.DataFrame(columns=[f"{disease}", "Genes"])
        with open(comm_path + f"results_100_{id}/seeds_{id}.txt", 'r') as file:
            i = 1
            for line in file:
                df._set_value(i, "Genes", line.rstrip())
                i += 1
        df.to_csv(path + f"output_tables/genes_comm_{id}.tsv", sep="\t", index=False)
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/genes_in_PA_communities.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()
    
list_genes_in_communities(orpha_file=orpha_names, list_comm=communities_100, list_ids_analyzed=list_ids_analyzed, comm_path=comm_path)
