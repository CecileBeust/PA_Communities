import networkx as nx
from scipy.stats import pearsonr
import pandas as pd
from pathlib import Path
import os
import argparse
import sys

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

data_folder = os.path.join(os.path.dirname(__file__), '..', '00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
cluster_output = os.path.join(data_folder, 'cluster_output_10_0.7.tsv')


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


def create_cluster_dico(cluster_file: str) -> dict :
    """Function to create a dictionary of disease
    communities clusters

    Args:
        cluster_file (str): name of the file containing
        the clusters assignments

    Returns:
        dict: the dictionary of disease communities
        clusters
    """
    df = pd.read_csv(cluster_file, sep="\t")
    dico_cluster_diseases = {}
    i = 0
    for cluster in df['cluster']:
        disease = df.iloc[i]['disease']
        if cluster not in dico_cluster_diseases.keys():
            dico_cluster_diseases[cluster] = [disease]
        else:
            dico_cluster_diseases[cluster] += [disease]
        i += 1
    return dico_cluster_diseases


dico_cluster_diseases = create_cluster_dico(cluster_output)
print(dico_cluster_diseases)

def filter_cluster(dico_cluster: dict) -> dict:
    """Function to filter a dictionary of 
    disease communities clusters to keep
    only clusters having at least 3 disease
    communities

    Args:
        dico_cluster (dict): the dico of
        disease commmunities clusters

    Returns:
        dict: the filtered dictionary 
    """
    filtered_dict = {}
    for cluster in dico_cluster:
        if len(dico_cluster[cluster]) >=3 :
            filtered_dict[cluster] = dico_cluster[cluster]
    return filtered_dict

filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(" ")
print(filtered_dico_cluster)

def extract_genes_clusters(filtered_dico_cluster: dict) -> None:
    """Function to extract the genes in a cluster of communities : generates
    an excel file gathering the gene composition of each cluster

    Args:
        filtered_dico_cluster (dict): dicitonary containing the assignement
        of communities in custers
    """
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/Pathways.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes.tsv", create_using = nx.Graph)
    for cluster in filtered_dico_cluster:
        print(cluster)
        diseases = []
        for disease in filtered_dico_cluster[cluster]:
            diseases.append(disease)
        genes_cluster = list()
        dico_genes_comm = dict()
        for disease in diseases:
            comm = comm_path + f"results_10_{disease}/seeds_{disease}.txt"
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
        for gene in genes_cluster:
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
        df = df.rename(columns={'Unnamed: 0': f"Cluster {cluster}"})
        df.to_csv(path + f"Clusters/genes_in_cluster_{cluster}.tsv", sep="\t")
        comm = df['Nb of communities in cluster']
        max_deg = df['max degree']
        sum_deg = df['sum degrees']
        corr_max, p_value_max = pearsonr(comm, max_deg)
        print('Correlation Nb comm / max deg:', corr_max)
        print('P-value:', p_value_max)
        corr_sum, p_value_sum = pearsonr(comm, sum_deg)
        print('Correlation Nb comm / sum deg:', corr_sum)
        print('P-value:', p_value_sum)
    os.mkdir(path + "Clusters")
    tsv_dir = Path(path + "Clusters")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "/genes_in_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()


extract_genes_clusters(filtered_dico_cluster)