import networkx as nx
from scipy.stats import pearsonr
import pandas as pd
from pathlib import Path
import os
import argparse
import sys

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../')
os.chdir(path)
print(path)

from utilities import create_cluster_dico, filter_cluster

data_folder = os.path.join(os.path.dirname(__file__), '..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')

# Argparse
parser = argparse.ArgumentParser(
    prog="analyse_genes_in_clusters.py", 
    description="functions to analyze the gene composition of a cluster of communities"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")

dico_cluster_diseases = create_cluster_dico(cluster_output)
print(f"Clusters: {dico_cluster_diseases}")

filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(f"Clusters containing at least 3 diseases: {filtered_dico_cluster}")

def extract_genes_clusters_and_degrees(filtered_dico_cluster: dict) -> None:
    """Function to extract the genes in a cluster of communities : generates
    an excel file gathering the gene composition of each cluster

    Args:
        filtered_dico_cluster (dict): dicitonary containing the assignement
        of communities in custers
    """
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    for cluster in filtered_dico_cluster:
        diseases = []
        for disease in filtered_dico_cluster[cluster]:
            diseases.append(disease)
        genes_cluster = list()
        dico_genes_comm = dict()
        for disease in diseases:
            comm = comm_path + f"results_100_{disease}/seeds_{disease}.txt"
            with open(comm, 'r') as file:
                for line in file:
                    gene = line.rstrip()
                    genes_cluster.append(gene)
                    if not gene in dico_genes_comm.keys():
                        dico_genes_comm[gene] = 1
                    else:
                        dico_genes_comm[gene] += 1
        i = 0
        df = pd.DataFrame(columns=["Gene", 'Nb of communities in cluster', 'degree PPI ntw', 'degree Pathways ntw', 'degree Coexp ntw', 'degree Complexes ntw', 'max degree', 'sum degrees'])
        for gene in set(genes_cluster):
            all_deg = []
            df._set_value(i, 'Gene', gene)
            df._set_value(i, 'Nb of communities in cluster', dico_genes_comm[gene])
            for elt in ppi.degree([gene]):
                df._set_value(i, 'degree PPI ntw', elt[1])
                all_deg.append(elt[1])
            for elt in pathways.degree([gene]):
                df._set_value(i, 'degree Pathways ntw', elt[1])
                all_deg.append(elt[1])
            for elt in coexp.degree([gene]):
                df._set_value(i, 'degree Coexp ntw', elt[1])
                all_deg.append(elt[1])
            for elt in complexes.degree([gene]):
                df._set_value(i, 'degree Complexes ntw', elt[1])
                all_deg.append(elt[1])
            df._set_value(i, 'max degree', max(all_deg))
            df._set_value(i, 'sum degrees', sum(all_deg))
            i += 1
        df = df.rename(columns={'Unnamed: 0': f"{cluster}"})
        df_sorted = df.sort_values(by=['Nb of communities in cluster'], ascending=False)
        df_sorted.to_csv(path + f"output_tables/genes_in_{cluster}_and_degrees.tsv", sep="\t", index=False)
        comm = df['Nb of communities in cluster']
        max_deg = df['max degree']
        sum_deg = df['sum degrees']
        corr_max, p_value_max = pearsonr(comm, max_deg)
        print('Correlation Nb comm / max deg:', corr_max)
        print('P-value:', p_value_max)
        corr_sum, p_value_sum = pearsonr(comm, sum_deg)
        print('Correlation Nb comm / sum deg:', corr_sum)
        print('P-value:', p_value_sum)
    tsv_dir = Path(path + "output_tables")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "output_tables/genes_in_clusters_and_degrees.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()


#extract_genes_clusters_and_degrees(filtered_dico_cluster)


def generate_supp_file(filtered_dico_cluster: dict) -> None:
    """Function to extract the genes in a cluster of communities : generates
    an excel file gathering the gene composition of each cluster

    Args:
        filtered_dico_cluster (dict): dicitonary containing the assignement
        of communities in custers
    """
    for cluster in filtered_dico_cluster:
        diseases = []
        for disease in filtered_dico_cluster[cluster]:
            diseases.append(disease)
        genes_cluster = list()
        for disease in diseases:
            comm = comm_path + f"results_100_{disease}/seeds_{disease}.txt"
            with open(comm, 'r') as file:
                for line in file:
                    gene = line.rstrip()
                    if not gene in genes_cluster:
                        genes_cluster.append(gene)
        i = 0
        df = pd.DataFrame(columns=[f"{cluster}", "Genes"])
        for gene in genes_cluster:
            df._set_value(i, 'Genes', gene)
            i += 1
        df.to_csv(path + f"output_tables/genes_in_{cluster}.tsv", sep="\t", index=False)
    tsv_dir = Path(path + "output_tables")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "output_tables/genes_in_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

generate_supp_file(filtered_dico_cluster=filtered_dico_cluster)