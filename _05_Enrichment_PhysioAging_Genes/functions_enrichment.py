import networkx as nx
import pandas as pd
import numpy as np
import scipy.stats as stats
import os

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

#############################################
###### LOAD AND PROCESS DATA ################
#############################################

def extract_seeds(file, list_ids_analyzed):
    orpha = pd.read_table(file)
    print(orpha)
    seeds = set()
    for index, row in orpha.iterrows():
        if str(row[0]) in list_ids_analyzed:
            seed = row[1].split(",")
            print(seed)
            for s in seed:
                seeds.add(s)
    return seeds

def load_geneage(file, seeds, nodes_ntw):
    # Load and filter Gene Age database
    df = pd.read_csv(file, sep=",")
    # remove genes without direct association with aging
    filtered_df = df.loc[(df['why'] != "upstream") & (df['why'] != "downstream") & (df['why'] != "putative") & (df['why'] != "functional") & (df['why'] != "downstream,putative") & (df['why'] != "functional,downstream") & (df['why'] != "functional,putative") & (df['why'] != "functional,upstream") & (df['why'] != "upstream,putative")]
    genes_aging = filtered_df['symbol'].to_list()
    print(f"{len(genes_aging)} aging genes")
    # remove seeds from aging genes to avoid overfitting
    genes_aging_wo_seeds = []
    for gene in genes_aging:
        if not gene in seeds:
            genes_aging_wo_seeds.append(gene)
    print(f"{len(genes_aging_wo_seeds)} aging genes after removing seeds genes")
    # check if aging genes are present in networks
    # remove aging gene not present in network
    genes_aging_wo_seeds_in_mx = list()
    for gene in genes_aging_wo_seeds:
        if gene in nodes_ntw:
            genes_aging_wo_seeds_in_mx.append(str(gene))
    print(f"{len(genes_aging_wo_seeds_in_mx)} aging genes without seeds in network")
    print(f"There are {len(set(genes_aging_wo_seeds_in_mx).intersection(nodes_ntw))} genage genes from the list which are in the multiplex network")
    return genes_aging_wo_seeds

def extract_genes_from_comm(comm_path: str, size: int, list_id_analyzed: list):
    dico_comm_nodes = {}
    for id in list_id_analyzed:
        community = comm_path + f"/results_{size}_{id}/seeds_{id}.txt"
        nodes_comm = []
        with open(community, 'r') as fi:
            for line in fi:
                gene = line.rsplit()
                nodes_comm.append(gene)
                dico_comm_nodes[id] = []
        for node in nodes_comm:
            dico_comm_nodes[id] += node
    return dico_comm_nodes

def extract_genes_from_cluster(comm_path: str, filtered_dico: dict, size: int):
    dico_clusters_nodes = dict()
    for cluster in filtered_dico:
        dico_clusters_nodes[cluster] = []
        for comm in filtered_dico[cluster]:
            community = comm_path + f"/results_{size}_{comm}/seeds_{comm}.txt"
            with open(community, 'r') as file:
                for line in file:
                    gene = line.rstrip()
                    if not gene in dico_clusters_nodes[cluster]:
                        dico_clusters_nodes[cluster] += [gene]
    return dico_clusters_nodes

def remove_seeds(gene_list: list, seeds_list: list) -> list:
    """Function which remove seeds nodes used to identify the communities
    from a list of genes

    Args:
        gene_list (list): list of genes to analyze
        seeds_list (list): list of seeds to remove

    Returns:
        list: new list of genes without seeds
    """
    new_list = []
    for gene in gene_list:
        if not gene in seeds_list:
            new_list.append(gene)
    return new_list


##################################################
######### MAPPING OF GENE IDENTIFIERS ############
##################################################

def create_mapping_file(rpkm_file: str) -> None:
    """Function to create a mapping file (Ensembl gene symbols to HUGO
    gene symbols) from the RPKM file of the study of Irizar et. al

    Args:
        rpkm_file (str): path to the RPKM file
    """
    df = pd.read_excel(rpkm_file)
    df = df.iloc[:, :2]
    print(df)
    df.to_csv(path + "/Data_PhysioAging/Mapping_Ensembl_GeneSymbol.txt", sep="\t", index=None)
    return path + "/Data_PhysioAging/Mapping_Ensembl_GeneSymbol.txt"

def getMappingDict(filePath: str, convertFrom: str, convertTo: str):
    df = pd.read_csv(filePath, sep="\t")
    df = df[[convertFrom, convertTo]]
    df = df.dropna(axis=0)
    if('ID' in df.columns):
        df['ID'] = df['ID'].astype(str)
        df['ID'] = df['ID'].astype(str)
    df = df.set_index(convertFrom)
    mappingDict = df.to_dict()[convertTo]
    return mappingDict

def map_ensembl_to_symbol(mapping_file_path: str, list_genes: list, From: str, To: str) -> list:
    """Function to map Ensembl genes identifiers to HUGO gene symbols

    Args:
        list_genes (list): genes to map
        From (str): identifiers to be mapped (Ensembl)
        To (str): identifiers to map (HUGO gene symbols)

    Returns:
        list: list of genes identified with HUGO genes symbols
    """
    mapping_dict = getMappingDict(
    filePath=mapping_file_path,
    convertFrom=From,
    convertTo=To
    )
    new_list = []
    for gene in list_genes:
        if gene in mapping_dict:
            new_list.append(mapping_dict[gene])
    return new_list

##################################################
######## ENRICHMENT FUNCTIONS ####################
##################################################

# Fischer's exact test : compare distributions of genes_comm and genes_aging_wo_seeds
def fisher(list1, list2, gene_pool):
    common_genes = set(list1).intersection(set(list2))
    print(f"{len(common_genes)} common genes")
    set1_unique = set(list1) - set(list2)
    print(f"{len(set1_unique)} genes in genes_comm but not in genes_aging_we_seeds")
    set2_unique = set(list2) - set(list1)
    print(f"{len(set2_unique)} genes in genes_aging_wo_seeds but not in genes_comm")

    not_set2 = gene_pool - len(set2_unique) - len(common_genes)
    #neither = gene_pool - len(common_genes) - len(set1_unique) - len(set2_unique)
    contingency_table = [
        [len(common_genes), len(set2_unique)],
        [len(set1_unique), not_set2]
    ]

    odds_ratio, p_value = stats.fisher_exact(contingency_table)

    print(f'Odds Ratio: {odds_ratio:.10f}, p-value: {p_value:.10f}')
    return p_value


def hypergeome(list1, list2, gene_pool):
    # Define the size of the gene pool (total number of genes)
    gene_pool_size = gene_pool

    # Define the number of genes in the two lists
    list1_size = len(list1)
    list2_size = len(list2)

    # Determine the number of genes that are common to both lists
    common_genes = set(list1).intersection(list2)
    common_genes_size = len(common_genes)

    # Define the number of genes to randomly sample from the gene pool
    sample_size = list1_size

    # Calculate the p-value using a hypergeometric test
    p_value = stats.hypergeom.sf(common_genes_size-1, gene_pool_size, list2_size, sample_size)

    print(f'p-value: {p_value:.4f}')
    return p_value
