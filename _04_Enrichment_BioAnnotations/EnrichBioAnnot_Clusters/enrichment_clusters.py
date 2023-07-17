import os
import argparse
from gprofiler import GProfiler
import ast
import pandas as pd
import sys
from pathlib import Path

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../..')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, build_communities_list, create_cluster_dico, filter_cluster, load_networks

data_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')
mapping_hugo_file = os.path.join(data_folder, 'mapping_HGNC_Ensembl.txt')

# Argparse
parser = argparse.ArgumentParser(
    prog="enrichment_clusters.py", 
    description="functions perform enrichment of clusters of gene communities"
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
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
# define background for hypergeometric test = number of unique nodes in the multiplex network
all_nodes_HUGO = list(load_networks(comm_path=comm_path))

def create_mapping_dict(mapping_file: str):
    mapping_dict = dict()
    df = pd.read_table(mapping_file, header=0)
    for index, row in df.iterrows():
        if str(row[0]) not in mapping_dict.keys():
            mapping_dict[str(row[0])] = str(row[1])
    return mapping_dict

mapping_dict = create_mapping_dict(mapping_file=mapping_hugo_file)

def map_all_nodes_mx(all_nodes: list, mapping_dict: dict):
    all_nodes_ensembl = list()
    for node in all_nodes:
        if node in mapping_dict:
            ensembl_node = mapping_dict[node]
            all_nodes_ensembl.append(ensembl_node)
        else:
            all_nodes_ensembl.append(node)
    return all_nodes_ensembl

all_nodes_ensembl = map_all_nodes_mx(all_nodes=all_nodes_HUGO, mapping_dict=mapping_dict)
print(len(all_nodes_ensembl))

pa_diseases = pd.read_csv(orpha_names, sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

dico_cluster_diseases = create_cluster_dico(cluster_output)
print(" ")
print(dico_cluster_diseases)

filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(" ")
print(f"Clusters containing at least 3 diseases: {filtered_dico_cluster}")

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
    print(genes)
    gp = GProfiler(user_agent='https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/gost/profile/', return_dataframe=True)
    enrich = gp.profile(organism='hsapiens', query=genes, sources=["GO:BP", "GO:CC", "REAC"], user_threshold=0.05, domain_scope='custom', significance_threshold_method='fdr', no_evidences=False, background=all_nodes_ensembl)
    print(enrich)
    enrich.to_csv(path + f"output_tables/enrich_bioannot_{cluster_id}.tsv", sep="\t")
    
def enrich_all_clusters(filtered_dico: dict) -> None:
    """Function to make the enrichment of all the clusters

    Args:
        filtered_dico (dict): teh dico of clusters and 
        the disease communities in it, containing only
        cluster having at least 3 disease communities
    """
    for cluster in filtered_dico.keys():
        print(cluster)
        nodes = select_nodes_from_cluster(comm_path, 100, filtered_dico, cluster)
        enrichment_cluster(cluster, nodes, 100)
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/enrich_bioannot_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

#enrich_all_clusters(filtered_dico_cluster)

def highlight_seeds(filtered_dico_cluster: dict, size: int, dico_code_disease: dict, dico_disease_seeds: dict) -> None:
    """Function which allows to add information to the enrichment files generated with enrich_all_clusters()
    It adds two columns to the enrichment tables : the seeds in the genes associated to each enrichment term, and
    the associated diseases

    Args:
        filtered_dico_cluster (dict): dictionary of clusters and their diseases
        size (int): number of iterations used for itRWR, reflecting the size of communities
        dico_code_disease (dict): dictionary containing PA diseases and their ORPHANET codes 
        dico_disease_seeds (dict): dictionary containing PA diseases and their associated genes
    """
    for cluster in filtered_dico_cluster.keys():
        print(" ")
        print("######################")
        print(cluster)
        seeds_cluster = list()
        diseases_cluster = list()
        for disease in filtered_dico_cluster[cluster]:
            diseases_cluster.append(str(disease))
            seed_disease = dico_disease_seeds[str(disease)]
            seeds_cluster.extend(seed_disease)
            seeds_cluster = list(set(seeds_cluster))
        print(diseases_cluster)
        enrichment_file = path + f"output_tables/enrich_bioannot_{cluster}.tsv"
        df = pd.read_csv(enrichment_file, sep="\t", header=0)
        df["seeds"] = 0
        df["diseases"] = 0
        i = 0
        for index, row in df.iterrows():
            intersection = ast.literal_eval(row['intersections'])
            seeds_genes = []
            diseases = []
            for gene in intersection:
                if gene in seeds_cluster:
                    seeds_genes.append(gene)
                    for key, value in dico_disease_seeds.items():
                        if gene in value:
                            if key in diseases_cluster:
                                diseases.append(key)
            if seeds_genes != []:
                df.iat[i, df.columns.get_loc('seeds')] = seeds_genes
            if diseases != []:
                dis_names = []
                for disease_code in diseases:
                    dis_name = dico_code_disease[disease_code]
                    dis_names.append(dis_name)
                df.iat[i, df.columns.get_loc('diseases')] = dis_names
            i += 1
        df = df.drop(columns=["evidences"])
        df.to_csv(path + f"output_tables/supp_files/enrich_bioannot_{cluster}.tsv", sep="\t")
    tsv_dir = Path(path + "output_tables/supp_files")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/supp_files/enrich_bioannot_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()
                
highlight_seeds(filtered_dico_cluster, 100, dico_code_disease, dico_disease_seeds)

def analyze_clusters(dico_clusters, dico_diseases_names):
    df = pd.DataFrame(columns=['Cluster', 'Diseases', 'Seeds'])
    print(df)
    i = 0
    for cluster in dico_clusters:
        print(cluster)
        diseases_codes = dico_clusters[cluster]
        diseases_names = []
        for code in diseases_codes:
            diseases_names.append(dico_diseases_names[str(code)])
        for disease, code in zip(diseases_names, diseases_codes):
            print(disease, code)
            seeds = dico_disease_seeds[str(code)]
            if cluster == 3:
                new_clust_id = 2
            elif cluster == 4:
                new_clust_id = 3
            elif cluster == 5:
                new_clust_id = 4
            elif cluster == 8:
                new_clust_id = 5
            elif cluster == 13:
                new_clust_id = 6
            else:
                new_clust_id = cluster
            new_row = [f"{new_clust_id}", disease, seeds]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
    print(df)
    df.to_csv(path + "output_tables/clusters_composition.tsv", sep="\t", index=False)