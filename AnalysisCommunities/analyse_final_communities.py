# Import modules
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

comm_path = "/home/cbeust/Landscape_PA/Github_Codes/PA_Communities/IDCommunity"

def get_list_orpha_names(path: str, pa_diseases: str, id_diseases_analysed: list) -> list:
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
    file = path + pa_diseases
    with open(file, 'r') as fi:
        for line in fi:
            disease_id = line.split("\t")[0].rstrip()
            disease_name = line.split("\t")[1].rstrip()
            for id in id_diseases_analysed:        
                if id == disease_id[6:]:
                    disease_name_list.append(disease_name)
    return disease_name_list

def create_dico_disease_seeds(path: str, orpha_seeds: str) -> tuple[dict, list]:
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
    file = path + orpha_seeds
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

(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "orpha_codes_PA.txt")
print(f"Dico diseases-seeds : {dico_disease_seeds}")
print(" ")
print(f"List diseases IDs : {list_id}")

def build_communities_list(path: str, list_id: list, size: int) -> tuple[list, list]:
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
        community = f"{path}/results_{size}_{id}/seeds_{id}.txt"
        # check is path exists = if the RWR has gave an output for this disease
        if os.path.exists(community):
            list_communities.append(community)
        # if not we store the id in a list of ids not analyzed
        else:
            not_analyzed.append(id)
    return (list_communities, not_analyzed)


(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)
print(f"Communities 10 : {communities_10}")
print(" ")
print(f"Diseases not analyzed (no seeds) : {not_analyzed}")

def jaccard_index(file1, file2) -> float:
    """Function to compute the Jaccard index
    between two communities stored in two
    files

    Args:
        file1 (str) : name of the first file
        file2 (str) : name of the second file

    Return:
        float : the value of the Jaccard index
        (intersection / union) between the 
        two file contents = the two
        communities
    """
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        contents1 = f1.read()
        contents2 = f2.read()
    
    words1 = contents1.split()
    words2 = contents2.split()
    
    set1 = set(words1)
    set2 = set(words2)
    intersection = len(set1.intersection(set2))
    union = len(set1) + len(set2) - intersection

    jaccard_index = intersection / union
    # semi jaccard index : need to set another argument
    # to tell on which set we want to compare
    # if nb == 1:
        #jaccard_index = intersection / len(set1)
    #elif nb == 2:
        #jaccard_index = intersection / len(set2)
    
    return jaccard_index


list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
print(f"There are {len(list_ids_analyzed)} diseases analyzed")
disease_name_list = get_list_orpha_names(path, "pa_orphanet_diseases.tsv", list_ids_analyzed)
print(" ")
print(f"Diseases names list : {disease_name_list}")

def build_similarity_matrix(nodes: list) -> tuple[np.ndarray, dict]:
    # Initialize similarity matrix
    similarity_matrix = np.zeros((len(nodes), len(nodes)))
    dico_jaccard = {}
    # Fill similarity matrix with Jaccard indexes
    for i in range(len(nodes)):
        for j in range(i, len(nodes)):
            similarity = jaccard_index(nodes[i], nodes[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
            if similarity not in dico_jaccard.keys():
                dico_jaccard[similarity] = 1
            else:
                dico_jaccard[similarity] += 1

    return (similarity_matrix, dico_jaccard)

(similarity_matrix_10, dico_jaccard) = build_similarity_matrix(communities_10)
print(" ")
print("Similarity matrix based on Jaccard index")
print(similarity_matrix_10)
print(dico_jaccard)

def plot_jaccard_distribution(dico_jaccard: dict, size: int) -> None:
    """Function to plot the distribution of
    the Jaccard indices among the communities

    Args:
        dico_jaccard (dict): dictionary with the values of the Jaccard indices
        size (int): number of iterations used to build the communities
    """
    mylist = [key for key, val in dico_jaccard.items() for _ in range(val)]
    plt.hist(mylist, bins=20, rwidth=0.9)
    plt.title(f"Distribution of Jaccard index for disease communities of size {size} + seeds")
    plt.show()

#plot_jaccard_distribution(dico_jaccard, 10)

def build_similarity_matrix_semi_jaccard(nodes: list):
    # Initialize similarity matrix
    similarity_matrix = np.zeros((len(nodes), len(nodes)))

    # Fill similarity matrix with Jaccard indexes
    for i in range(len(nodes)):
        for j in range(i, len(nodes)):
            similarity_1 = jaccard_index(nodes[i], nodes[j], 1)
            similarity_matrix[i, j] = similarity_1
            similarity_2 = jaccard_index(nodes[i], nodes[j], 2)
            similarity_matrix[j, i] = similarity_2
    print(similarity_matrix)
    return similarity_matrix

#sim_matrix_semi_jaccard_50 = build_similarity_matrix_semi_jaccard(communities_50)

def build_sim_matrix_seeds(dico_seeds, list_id):
    similarity_matrix = np.zeros((len(list_id), len(list_id)))
    for i in range(len(list_id)):
        id_1 = list_id[i]
        for j in range(i, len(list_id)):
            id_2 = list_id[j]
            seeds_1 = dico_seeds[id_1]
            seeds_2 = dico_seeds[id_2]
            similarity = jaccard_index(seeds_1, seeds_2)
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    return similarity_matrix
