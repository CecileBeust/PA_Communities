import os
import argparse
from gprofiler import GProfiler
import ast
import pandas as pd
import sys
from pathlib import Path
import numpy as np

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
mapping_hugo_file = os.path.join(data_folder, 'mapping_HGNC_Ensembl_170723.txt')
reac_genes = "../Build_background/REAC_genes.csv"
gobp_genes = "../Build_background/GOBP_genes.csv"
gocc_genes = "../Build_background/GOCC_genes.csv"

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
all_nodes_HUGO = list(load_networks(comm_path=comm_path))

def compute_background_annot_mx(genes_file_annot: str, multiplex_nodes: list) -> list():
    df = pd.read_csv(genes_file_annot, sep=",", header=0)
    genes_annot = df["Gene"].to_list()
    background = set(genes_annot).intersection(set(multiplex_nodes))
    return list(background)

background_gobp = compute_background_annot_mx(genes_file_annot=gobp_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg go bp hugo: {len(background_gobp)}")
background_gocc = compute_background_annot_mx(genes_file_annot=gocc_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg go cc hugo: {len(background_gocc)}")
background_reac = compute_background_annot_mx(genes_file_annot=reac_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg reac hugo: {len(background_reac)}")

pa_diseases = pd.read_csv(orpha_names, sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(" ")
print(f"Dico diseases code: {dico_code_disease}")

dico_cluster_diseases = create_cluster_dico(cluster_output)
print(" ")
print(f"Dico clusters diseases: {dico_cluster_diseases}")

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

def enrichment_cluster(source: str, cluster_id: int, gene_set: set, background: list) -> None:
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
    gp = GProfiler(
        user_agent='https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/gost/profile/', 
        return_dataframe=True
        )
    enrich = gp.profile(
        organism='hsapiens', 
        query=genes,
        sources=[source],
        domain_scope='custom',
        user_threshold=0.05,
        significance_threshold_method='false_discovery_rate', 
        background=background,
        no_evidences=False,
        )
    source = source.replace(":", "")
    enrich.to_csv(path + f"output_tables/enrich_{source}_{cluster_id}.tsv", sep="\t")

def enrich_all_clusters(filtered_dico: dict) -> None:
    """Function to make the enrichment of all the clusters

    Args:
        filtered_dico (dict): teh dico of clusters and 
        the disease communities in it, containing only
        cluster having at least 3 disease communities
    """
    for cluster in filtered_dico.keys():
        nodes = select_nodes_from_cluster(comm_path=comm_path, size=100, dico_cluster=filtered_dico, cluster_id=cluster)
        enrichment_cluster(source="GO:BP", cluster_id=cluster, gene_set=nodes, background=background_gobp)
        enrichment_cluster(source="GO:CC", cluster_id=cluster, gene_set=nodes, background=background_gocc)
        enrichment_cluster(source="REAC", cluster_id=cluster, gene_set=nodes, background=background_reac)
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/enrich_bioannot_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

enrich_all_clusters(filtered_dico=filtered_dico_cluster)

def merge_enrichment_files(filetered_dico_cluster: dict, sources: list) -> None:
    for cluster in filetered_dico_cluster:
        gobp = pd.read_csv(f"output_tables/enrich_{sources[0]}_{cluster}.tsv", sep="\t", header=0)
        gocc = pd.read_csv(f"output_tables/enrich_{sources[1]}_{cluster}.tsv", sep="\t", header=0)
        reac = pd.read_csv(f"output_tables/enrich_{sources[2]}_{cluster}.tsv", sep="\t", header=0)
        all_annot = pd.concat([gobp, gocc, reac], axis=0)
        assert len(all_annot.index) == len(gobp.index) + len(gocc.index) + len(reac.index)
        all_annot = all_annot.drop(all_annot.columns[0], axis=1)
        print(all_annot)
        all_annot_sorted = all_annot.sort_values(by=["p_value"], ascending=True)
        all_annot_sorted.rename(columns={'p_value':'Corrected p_value'}, inplace=True)
        all_annot_sorted.to_csv(f"output_tables/enrich_bioannot_{cluster}.tsv", sep="\t", header=True, index=False)
        os.remove(f"output_tables/enrich_{sources[0]}_{cluster}.tsv")
        os.remove(f"output_tables/enrich_{sources[1]}_{cluster}.tsv")
        os.remove(f"output_tables/enrich_{sources[2]}_{cluster}.tsv")
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/enrich_bioannot_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

merge_enrichment_files(filetered_dico_cluster=filtered_dico_cluster, sources=["GOBP", "GOCC", "REAC"])

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
        print(cluster)
        seeds_cluster = list()
        diseases_cluster = list()
        for disease in filtered_dico_cluster[cluster]:
            diseases_cluster.append(str(disease))
            seed_disease = dico_disease_seeds[str(disease)]
            seeds_cluster.extend(seed_disease)
            seeds_cluster = list(set(seeds_cluster))
        print(diseases_cluster)
        print(seeds_cluster)
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
        df.to_csv(path + f"output_tables/supp_files/enrich_bioannot_{cluster}.tsv", sep="\t", header=True, index=False)
    tsv_dir = Path(path + "output_tables/supp_files/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + f"output_tables/supp_files/enrich_bioannot_clusters.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()
                
highlight_seeds(filtered_dico_cluster, 100, dico_code_disease, dico_disease_seeds)