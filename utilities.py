# Import modules
import os
import pandas as pd
import networkx as nx
import numpy as np

def get_list_orpha_names(pa_diseases: str, id_diseases_analysed: list) -> list:
    """Function to get a list of PA diseases names
    from a file of diseases + their IDs and a list 
    of identifiers to analyze

    Args:
        path (str) : path of the working directory
        pa_diseases (str): _description_
        id_diseases_analysed (list): list of disease
        identifiers that we want to analyze

    Returns:
        list: list of disease names to analyze
    """
    disease_name_list = []
    file = pa_diseases
    with open(file, 'r') as fi:
        for line in fi:
            disease_id = line.split("\t")[0].rstrip()
            disease_name = line.split("\t")[1].rstrip()
            for id in id_diseases_analysed:        
                if id == disease_id[6:]:
                    disease_name_list.append(disease_name)
    return disease_name_list

def create_dico_disease_seeds(orpha_codes: str) -> tuple[dict, list]:
    """Function to create a dictionary of PA
    diseases and their associated causative
    genes (seeds)

    Args:
        path (str) : path of the working directory
        orpha_seeds (str): name of the file 
        containing the diseases IDs and their
        associated causative genes

    Returns:
        dict: dictionnary containing disease
        identifiers as keys and their seeds as
        values
        list : the list of disease identifiers
    """
    dico_seeds = {}
    list_disease = []
    file = orpha_codes
    with open(file, 'r') as fi:
        for line in fi:
            # separate diseases from seeds in the input file
            disease = line.split("\t")[0]
            list_disease.append(disease)
            # split seeds
            seeds_rsplit = line.split("\t")[1].rsplit()
            seeds = [genes.split(",") for genes in seeds_rsplit]
            # initialize key in dico for disease
            dico_seeds[disease] = []
            # writing one seeds file for each set of seeds
            # we take the orpha code of the disease to name the seeds files
            for list_genes in seeds:
                for genes in list_genes:
                    # add set of seeds in dico
                    dico_seeds[disease].append(genes)
        return (dico_seeds, list_disease)
    
def build_communities_list(comm_path: str, list_id: list, size: int) -> tuple[list, list]:
    """Function which builds a list of communities by gathering
    all communities stored in different forlders identified by
    identifiers after itRWR

    Args:
        path (str): the path to the directory where RWR output
        folder containing communities are stored
        list_id (list): list of ids for the communities
        (corresponding to the ORPHANET codes of the diseases)
        size (int) : the number of iterations used to build
        the communities

    Returns:
        list: the list of communities names
        list: list of IDs of diseases that we do not analyze
        (because they have no seeds)
    """
    list_communities = []
    not_analyzed = []
    for id in list_id:
        # set path for community corresponding to the id
        community = comm_path + f"results_{size}_{id}/seeds_{id}.txt"
        # check is path exists = if the RWR has gave an output for this disease
        if os.path.exists(community):
            list_communities.append(community)
        # if not we store the id in a list of ids not analyzed
        else:
            not_analyzed.append(id)
    return (list_communities, not_analyzed)

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
    print(" ")
    print(filtered_dict)
    # Rename clusters by order of appearance on the clustermap
    filtered_dict["cluster_1"] = filtered_dict.pop(1)
    filtered_dict["cluster_2"] = filtered_dict.pop(3)
    filtered_dict["cluster_3"] = filtered_dict.pop(4)
    filtered_dict["cluster_4"] = filtered_dict.pop(5)
    filtered_dict["cluster_5"] = filtered_dict.pop(8)
    filtered_dict["cluster_6"] = filtered_dict.pop(13)
    print(" ")
    print(filtered_dict)
    return filtered_dict

def load_networks(comm_path: str):
    # load networks
    ppi = nx.read_edgelist(comm_path + "/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "/multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)

    # get noeds in networks
    ppi_nodes = ppi.nodes()
    pat_nodes = pathways.nodes()
    coexp_nodes = coexp.nodes()
    complexes_nodes = complexes.nodes()
    #all_nodes = list(set(list(ppi_nodes) + list(pat_nodes) + list(coexp_nodes) + list(complexes_nodes) + list(diseases_nodes)))
    all_nodes = np.unique(list(ppi_nodes) + list(pat_nodes) + list(coexp_nodes) + list(complexes_nodes))
    print(f"{len(all_nodes)} nodes in the multiplex network")
    return all_nodes