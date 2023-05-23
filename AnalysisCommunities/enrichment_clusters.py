import os
import argparse
from enrichment_communities import create_dico_disease_seeds, build_communities_list
from gprofiler import GProfiler
import ast
import pandas as pd

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# Argparse
parser = argparse.ArgumentParser(
    prog="cluster_communities.py", 
    description="functions cluster communities based on Jaccard index"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")


(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "data/orpha_codes_PA.txt")
(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)

pa_diseases = pd.read_csv(path + 'data/pa_orphanet_diseases.tsv', sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

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

dico_cluster_diseases = create_cluster_dico("data/cluster_output_10_0.7.tsv")
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
print(filtered_dico_cluster)

def select_seeds_from_cluster(dico_disease_seeds: dict, dico_cluster: dict, cluster_id: int) -> set:
    """Function to select the seeds genes from a cluster

    Args:
        dico_disease_seeds (dict) : dico containing the diseases
        and their associated genes (seeds)
        dico_cluster (dict): the dico of disease communities
        clusters
        cluster_id (int): the ID (number) of the cluster we want
        to select the seeds from

    Returns:
        set: set of seeds from the cluster
    """
    seeds = set()
    for disease in dico_cluster[cluster_id]:
        seeds_disease = dico_disease_seeds[str(disease)]
        for seed in seeds_disease:
            seeds.add(seed)
    return seeds


def select_nodes_from_cluster(comm_path: str, size: int, dico_cluster: dict, cluster_id: int) -> set:
    """Function to select the nodes from a disease community cluster

    Args:
        comm_path (str) : working directory containing the results of the
        community identification
        size (int): the number of iterations chosen to build the communities
        dico_cluster (dict): dico containing the clusters and the diseases communities in it
        cluster_id (int): the cluster we want to select the nodes from

    Returns:
        set: set of nodes of the cluster
    """
    nodes = set()
    for disease in dico_cluster[cluster_id]:
        with open(comm_path + f"results_{size}_{str(disease)}/seeds_{str(disease)}.txt", 'r') as community:
            for line in community:
                nodes.add(line.rstrip())
    return nodes


def enrichment_cluster(cluster_id: int, gene_set: set, size: int):
    """Function to enrich a disease communities cluster wth g:Profiler

    Args:
        cluster_id (int) : the cluster to enrich
        gene_set (set) : set of genes in the cluster
        size (int) : the number of iterations chosen to 
        build the communities
    Return:
        None
    """
    genes = list(gene_set)
    gp = GProfiler(return_dataframe=True)
    enrich = gp.profile(organism='hsapiens', query=genes, no_evidences=False)
    print(enrich)
    enrich.to_csv(path + f"cluster_{size}_{cluster_id}.tsv", sep="\t")
    
def enrich_all_clusters(filtered_dico: dict) -> None:
    """Function to make the enrichment of all the clusters

    Args:
        filtered_dico (dict): teh dico of clusters and 
        the disease communities in it, containing only
        cluster having at least 3 disease communities
    """
    for cluster in filtered_dico.keys():
        print(cluster)
        nodes = select_nodes_from_cluster(100, filtered_dico, cluster)
        enrichment_cluster(cluster, nodes, 100)

#enrich_all_clusters(filtered_dico_cluster_10_0_7)

def highlight_seeds(filtered_dico_cluster: dict, size: int, dico_code_disease: dict, dico_disease_seeds: dict):
    seeds = set()
    for list_gene in dico_disease_seeds.values():
        for gene in list_gene:
            seeds.add(gene)
    for cluster in filtered_dico_cluster.keys():
        enrichment_file = path + f"cluster_{size}_{cluster}.tsv"
        df = pd.read_csv(enrichment_file, sep="\t")
        df["seed"] = 0
        df["disease"] = 0
        i = 0
        for index, row in df.iterrows():
            intersection = ast.literal_eval(row['intersections'])
            seeds_genes = []
            diseases = []
            for gene in intersection:
                if gene in seeds:
                    seeds_genes.append(gene)
                    for key, value in dico_disease_seeds.items():
                        if gene in value:
                            diseases.append(key)
            if seeds_genes != []:
                df._set_value(i, 'seed', str(seeds_genes))
            if diseases != []:
                dis_names = []
                for disease_code in diseases:
                    dis_name = dico_code_disease[disease_code]
                    dis_names.append(dis_name)
                df._set_value(i, 'disease', str(dis_names)) 
            i += 1
        df.to_csv(path + f"cluster_{size}_{cluster}_with_genes.tsv", sep="\t")
                
#highlight_seeds(filtered_dico_cluster, 10, dico_code_disease, dico_disease_seeds)