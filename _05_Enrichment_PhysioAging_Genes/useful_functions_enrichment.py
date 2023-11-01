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
    orpha = pd.read_table(file, header=None)
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
    return genes_aging_wo_seeds_in_mx

def load_geneage_keep_seeds(file, seeds, nodes_ntw):
    # Load and filter Gene Age database
    df = pd.read_csv(file, sep=",")
    # remove genes without direct association with aging
    filtered_df = df.loc[(df['why'] != "upstream") & (df['why'] != "downstream") & (df['why'] != "putative") & (df['why'] != "functional") & (df['why'] != "downstream,putative") & (df['why'] != "functional,downstream") & (df['why'] != "functional,putative") & (df['why'] != "functional,upstream") & (df['why'] != "upstream,putative")]
    genes_aging = filtered_df['symbol'].to_list()
    print(f"{len(genes_aging)} aging genes")
    # check if aging genes are present in networks
    # remove aging gene not present in network
    genes_aging_in_mx = list()
    for gene in genes_aging:
        if gene in nodes_ntw:
            genes_aging_in_mx.append(str(gene))
    print(f"{len(genes_aging_in_mx)} aging genes without seeds in network")
    print(f"There are {len(set(genes_aging_in_mx).intersection(nodes_ntw))} genage genes from the list which are in the multiplex network")
    return genes_aging_in_mx

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

def extract_list_genes_DEG_physio_aging(file: str) -> tuple[list, list, list]:
    """Function to extract the lists of genes up, down or up and down-regulated
    during physiological aging described in the study of Irizar et. al

    Args:
        enrichment_file (str): path to the file containing the list of genes

    Returns:
        tuple[list, list, list]: list of genes up-regulated during physiological 
        aging in the studied tissue, list of genes down-regulated, and list of 
        genes both up and down-regulated. 
    """
    up_genes = []
    down_genes = []
    other_genes = []
    is_up = None
    with open(file, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            line = line.strip()
            if line == "Up:":
                is_up = True
            elif line == "Down:":
                is_up = False
            elif line == "Other:":
                 is_up = None
            elif is_up is True:
                up_genes.append(line)
            elif is_up is False:
                down_genes.append(line)
            else:
                other_genes.append(line)
    print("Up-regulated genes:", len(up_genes))
    print("Down-regulated genes:", len(down_genes))
    print("Other genes :", len(other_genes))
    return up_genes, down_genes, other_genes

def create_filtered_enrichment_lists_physio_aging_DEG(file: str, mapping_file_path: str, seeds_list: list, is_up: int, all_nodes: list) -> list:
    """Function to create the lists of genes up, down or both up and down regulated
    during physiological aging for the enrichments analysis. The seeds are removed from
    the list. 

    Args:
        file (str): path to file containing the list of genes
        seeds_list (list): list of seeds nodes
        is_up (int): 1 if we perform enrichment with genes up-regulated, 0 if we 
        performe enrichment with genes down-regulated
        all_nodes (list): list of all nodes in the multiplex network

    Returns:
        list: list of genes for the enrichment analysis
    """
    # extract lists of genes from file
    genes_up, genes_down, other_genes = extract_list_genes_DEG_physio_aging(file)

    # If we want to analyze UP regulated genes
    if is_up == 1:
        # map to gene symbol
        genes_up_GS = map_ensembl_to_symbol(mapping_file_path, genes_up, 'ID', 'external_gene_id')
        print(f"{len(genes_up_GS)} genes UP")
        # remove seeds
        gene_up_wo_seeds = remove_seeds(genes_up_GS, seeds_list)
        print(f"{len(gene_up_wo_seeds)} genes UP without the seeds")
        # check if aging genes are present in networks
        genes_up_wo_seeds_in_ntw = []
        for gene in gene_up_wo_seeds:
            if gene in all_nodes:
                genes_up_wo_seeds_in_ntw.append(gene)
        print(f"{len(genes_up_wo_seeds_in_ntw)} genes UP present in multiplex (WO seeds)")
        return genes_up_wo_seeds_in_ntw
    
    # If we want to analyze DOWN regulated genes
    elif is_up == 0:
        # map to gene names
        genes_down_GS = map_ensembl_to_symbol(mapping_file_path, genes_down, 'ID', 'external_gene_id')
        print(f"{len(genes_down_GS)} genes DOWN")
        # remove seeds
        gene_down_wo_seeds = remove_seeds(genes_down_GS, seeds_list)
        print(f"{len(gene_down_wo_seeds)} genes DOWN without seeds")
        # check if aging genes are present in networks
        genes_down_wo_seeds_in_ntw = []
        for gene in gene_down_wo_seeds:
            if gene in all_nodes:
                genes_down_wo_seeds_in_ntw.append(gene)
        print(f"{len(genes_down_wo_seeds_in_ntw)} genes DOWN present in multiplex (WO seeds)")
        return genes_down_wo_seeds_in_ntw
    

def create_filtered_enrichment_lists_physio_aging_DEG_keep_seeds(file: str, mapping_file_path: str, seeds_list: list, is_up: int, all_nodes: list) -> list:
    """Function to create the lists of genes up, down or both up and down regulated
    during physiological aging for the enrichments analysis. The seeds are removed from
    the list. 

    Args:
        file (str): path to file containing the list of genes
        seeds_list (list): list of seeds nodes
        is_up (int): 1 if we perform enrichment with genes up-regulated, 0 if we 
        performe enrichment with genes down-regulated
        all_nodes (list): list of all nodes in the multiplex network

    Returns:
        list: list of genes for the enrichment analysis
    """
    # extract lists of genes from file
    genes_up, genes_down, other_genes = extract_list_genes_DEG_physio_aging(file)

    # If we want to analyze UP regulated genes
    if is_up == 1:
        # map to gene symbol
        genes_up_GS = map_ensembl_to_symbol(mapping_file_path, genes_up, 'ID', 'external_gene_id')
        print(f"{len(genes_up_GS)} genes UP")
        # check if aging genes are present in networks
        genes_up_in_ntw = []
        for gene in genes_up_GS:
            if gene in all_nodes:
                genes_up_in_ntw.append(gene)
        print(f"{len(genes_up_in_ntw)} genes UP present in multiplex (WO seeds)")
        return genes_up_in_ntw
    
    # If we want to analyze DOWN regulated genes
    elif is_up == 0:
        # map to gene names
        genes_down_GS = map_ensembl_to_symbol(mapping_file_path, genes_down, 'ID', 'external_gene_id')
        print(f"{len(genes_down_GS)} genes DOWN")
        # check if aging genes are present in networks
        genes_down_in_ntw = []
        for gene in genes_down_GS:
            if gene in all_nodes:
                genes_down_in_ntw.append(gene)
        print(f"{len(genes_down_in_ntw)} genes DOWN present in multiplex (WO seeds)")
        return genes_down_in_ntw

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
def fisher(list1: list, list2: list, gene_pool: list):
    common_genes = set(list1).intersection(set(list2))
    set1_unique = set(list1) - set(list2)
    set2_unique = set(list2) - set(list1)

    not_set2 = len(gene_pool) - len(set2_unique) - len(common_genes)
    #neither = gene_pool - len(common_genes) - len(set1_unique) - len(set2_unique)
    contingency_table = [
        [len(common_genes), len(set2_unique)],
        [len(set1_unique), not_set2]
    ]

    odds_ratio, p_value = stats.fisher_exact(contingency_table)

    return p_value


def hypergeome(list1: list, list2: list, gene_pool: list):
    # Define the size of the background
    background = len(gene_pool)

    # Define the number of genes in the two lists
    list1_size = len(list1)
    list2_size = len(list2)

    # Determine the number of genes that are common to both lists
    common_genes = set(list1).intersection(set(list2))
    common_genes_size = len(common_genes)

    # Calculate the p-value using a hypergeometric test
    p_value = stats.hypergeom.sf(common_genes_size-1, background, list2_size, list1_size)

    return p_value