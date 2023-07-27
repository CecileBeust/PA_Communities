"""
Functions to perform a clustering of 
gene communities based on Jaccard index
"""

# Import modules and define path variable
import os
import sys
import shutil
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import argparse

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, get_list_orpha_names, build_communities_list, create_cluster_dico, filter_cluster

data_folder = os.path.join(os.path.dirname(__file__), '..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')

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

pa_diseases = pd.read_csv(orpha_names, sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

# Define the dico of diseases and their seeds + list of ORPHANET identifiers\
(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)

# Build the list of communities
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
print(f"Number of communities : {len(communities_100)}")
print(" ")
print(f"Diseases not analyzed : {not_analyzed}")
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
print(f"There are {len(list_ids_analyzed)} diseases analyzed")
disease_name_list = get_list_orpha_names(orpha_names, list_ids_analyzed)
print(" ")
print(f"Diseases names list : {disease_name_list}")

# build a list of disease names
df_disease = pd.read_csv(orpha_names, sep="\t", header=None)
disease_names = []
for id in list_ids_analyzed:
    for row in df_disease.iterrows():
        if f"ORPHA:{id}" == row[df_disease.columns[1]][0]:
            disease_names.append(str(row[df_disease.columns[1]][1]))

# Define a list of shorten disease names for better further visualization
shorter_disease_names = ['FG syndrome type 1', 
                         'Ablepharon macrostomia syndrome', 
                         'Cerebrotendinous xanthomatosis', 
                         'Williams syndrome', 
                         'Autosomal dominant cutis laxa', 
                         'Cockayne syndrome type 3', 
                         'Cockayne syndrome type 2', 
                         'Cockayne syndrome type 1', 
                         'Werner syndrome', 
                         'Mandibuloacral dysplasia with type B lipodystrophy', 
                         'Mandibuloacral dysplasia with type A lipodystrophy', 
                         'Waardenburg-Shah syndrome', 
                         'Waardenburg syndrome type 2', 
                         'Waardenburg syndrome type 1', 
                         'Seckel syndrome', 
                         'Griscelli syndrome type 2', 
                         'Griscelli syndrome type 1', 
                         'Atypical Werner syndrome', 
                         'COG7-CDG', 
                         'ALG8-CDG', 
                         'Acquired partial lipodystrophy', 
                         'Rabson-Mendenhall syndrome', 
                         'Pseudoxanthoma elasticum', 
                         'B4GALT7-related spondylodysplastic Ehlers-Danlos syndrome', 
                         'Hutchinson-Gilford progeria syndrome', 
                         'Laron syndrome', 
                         'Noonan syndrome with multiple lentigines', 
                         'Aquagenic palmoplantar keratoderma', 
                         'Pierpont syndrome', 
                         'KID syndrome', 
                         'Keppen-Lubinsky syndrome', 
                         'Autosomal recessive cerebellar ataxia due to STUB1 deficiency', 
                         'LMNA-related cardiocutaneous progeria syndrome', 
                         'Autosomal recessive cutis laxa type 2, classic type', 
                         'Wiedemann-Rautenstrauch syndrome', 
                         'Vogt-Koyanagi-Harada disease', 
                         'Neuroectodermal melanolysosomal disease', 
                         'Trichothiodystrophy', 
                         'Hoyeraal-Hreidarsson syndrome', 
                         'SHORT syndrome', 
                         'Nicolaides-Baraitser syndrome', 
                         'Progeroid syndrome, Petty type', 
                         'Classical Ehlers-Danlos syndrome', 
                         'Vascular Ehlers-Danlos syndrome', 
                         'Wrinkly skin syndrome', 
                         'MASSFDHH syndrome', 
                         'Autosomal semi-dominant severe lipodystrophic laminopathy', 
                         'Lenz-Majewski hyperostotic dwarfism', 
                         'Acral peeling skin syndrome', 
                         'COG5-CDG', 
                         'Acrogeria', 
                         'Xeroderma pigmentosum-Cockayne syndrome complex', 
                         'Geroderma osteodysplastica', 
                         'GAPO syndrome', 
                         'Myoclonic-astatic epilepsy', 
                         'Dermatosparaxis Ehlers-Danlos syndrome', 
                         'Thanatophoric dysplasia type 1', 
                         'Focal facial dermal dysplasia type III', 
                         'Dyskeratosis congenita', 
                         'PDNCDLWSH disease', 
                         'Cataract-intellectual disability-hypogonadism syndrome', 
                         'SOLAMEN Syndrome', 
                         'Cardiofaciocutaneous syndrome', 
                         'Branchioskeletogenital syndrome', 
                         'Branchio-oculo-facial syndrome', 
                         'Transaldolase deficiency', 
                         'Ataxia-telangiectasia'
                         ]

# create a dictionary of diseases ORPAHNET codes and their shorter names
assert len(list_ids_analyzed) == len(shorter_disease_names)
dico_id_shorter_names = dict()
for code, name in zip(list_ids_analyzed, shorter_disease_names):
    dico_id_shorter_names[code] = name
print(dico_id_shorter_names)

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
    
    return jaccard_index

def build_similarity_matrix(comm: list) -> tuple[np.ndarray, dict]:
    """Function allowing to build a similarity matrix of
    communities based on Jaccard index

    Args:
        comm (list): list of communities 

    Returns:
        tuple[np.ndarray, dict]: the similarity matrix and
        a dictionary of freauency of Jaccard indices values
    """
    # Initialize similarity matrix
    similarity_matrix = np.zeros((len(comm), len(comm)))
    dico_jaccard = {}
    # Fill similarity matrix with Jaccard indexes
    for i in range(len(comm)):
        for j in range(i, len(comm)):
            similarity = jaccard_index(comm[i], comm[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
            if similarity not in dico_jaccard.keys():
                dico_jaccard[similarity] = 1
            else:
                dico_jaccard[similarity] += 1

    return (similarity_matrix, dico_jaccard)

# Build similarity matrix of communities based on Jaccard index
(similarity_matrix_100, dico_jaccard) = build_similarity_matrix(comm=communities_100)
print(" ")
print("Similarity matrix based on Jaccard index")
print(similarity_matrix_100)

def build_distance_matrix(sim_matrix) -> np.ndarray:
    """Function to build a distance 
    matrix from a similarity matrix

    Args:
        sim_matrix (np.ndarray): a similaruty
        matrix

    Returns:
        np.ndarray: a distance matrix
    """
    new_matrix = np.zeros((len(sim_matrix), len(sim_matrix)))
    for i in range(len(sim_matrix)):
        for j in range(i, len(sim_matrix)):
            val = sim_matrix[i][j]
            new_matrix[i][j] = 1-val
            new_matrix[j][i] = 1-val
    return new_matrix

# Build a distance matrix from the similarity matrix
distance_matrix_100 = build_distance_matrix(sim_matrix=similarity_matrix_100)

def cluster_dendrogram(matrix: np.ndarray, cutoff: float, size: int, dico_id_shorter_names: dict):
    """
    Function to plot a hierarchichal clustering dendrogram from a distance matrix

    Args:
        matrix (np.ndarray) : a distance matrix
        cutoff (float) : a value to cut the tree
        size (int) : the number of iterations chose
        to build the communities
        dic_shorter_names (dict) : dictionary containing diseases
        ORPHANET codes as keys and their names (shortened) as values
    Return:
        None
    """
    # transform distance matrix to flat array, only keep the upper part of the matrix
    array = matrix[np.triu_indices(len(matrix), k=1)]
    # linkage matrix 
    Z = linkage(array, 'average')
    
    # plot dendrogram
    fig, ax = plt.subplots(figsize=(15, 15))
    plt.title("Hierarchichal clustering of PA diseases-associated communities")
    dn = dendrogram(Z, labels=list(dico_id_shorter_names.values()), leaf_font_size=7, ax=ax, color_threshold=None, orientation='left')
    plt.savefig(path + "output_figures/Dendrogram_PA_diseases.png", bbox_inches='tight')
    # define clusters assignments
    assignments = fcluster(Z, t=cutoff, criterion='distance')
    print(f"Cluster assignements : {assignments}") 
    print(f"Number of clusters : {max(assignments)}")
    print(len(list_ids_analyzed))
    
    # export cluster assignment to dataframe
    cluster_output = pd.DataFrame({'disease': list_ids_analyzed, 'cluster':assignments})
    # export cluster output to file
    cluster_output.to_csv(f"cluster_output_{size}_{cutoff}.tsv", sep="\t")
    source_file = f"cluster_output_{size}_{cutoff}.tsv"
    destination_folder = path + "../_00_data/"
    shutil.copy(source_file, destination_folder)
    
    # add cutoff line to dendrogram
    plt.axvline(x=cutoff, c='r')
    plt.tight_layout()
    
    # export distance matrix to dataframe for clustermap
    df = pd.DataFrame(matrix, columns=list(dico_id_shorter_names.values()), index=list(dico_id_shorter_names.values()))
    print("Dataframe of distance matrix")
    print(df)
    
    sns.set(font_scale=0.55)
    sns.set_style("white")
    cluster_colors = sns.color_palette("hls", n_colors=len(set(cluster_output['cluster'])))
    cluster_color_dict = dict(zip(sorted(set(cluster_output['cluster'])), cluster_colors))
    cluster_sizes = cluster_output.groupby('cluster').size()
    row_colors = [cluster_color_dict[label] if cluster_sizes[label] >= 3 else 'gray' for label in cluster_output['cluster']]
    clust = sns.clustermap(
        df, 
        annot=False,
        metric="euclidean",
        method="average",
        cmap="magma",
        row_linkage=Z, 
        col_linkage=Z,
        figsize=(10,10),
        xticklabels = True,
        yticklabels = False,
        square=True,
        linewidths=2,
        col_colors=row_colors
                )
    # mask lower part of diagonal
    mask = np.tril(np.ones_like(matrix), k=-1)
    values = clust.ax_heatmap.collections[0].get_array().reshape(matrix.shape)
    new_values = np.ma.array(values, mask=mask)
    clust.ax_heatmap.collections[0].set_array(new_values)
    clust.ax_row_dendrogram.set_visible(False)
    plt.title("Clustering of Premature Aging diseases based on the Jaccard index similarity of their identified communities", loc="left")
    plt.savefig(path + '/output_figures/clustering_PA_diseases.png')
    plt.show()

cluster_dendrogram(matrix=distance_matrix_100, cutoff=0.7, size=100, dico_id_shorter_names=dico_id_shorter_names)
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')

dico_cluster_diseases = create_cluster_dico(cluster_output)
print(" ")
print(f"Clusters: {dico_cluster_diseases}")

filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(" ")
print(f"Clusters containing at least 3 diseases: {filtered_dico_cluster}")

def analyze_clusters(dico_disease_seeds: dict, filtered_dico_cluster: dict, dico_code_disease: dict) -> None:
    df = pd.DataFrame(columns=['Cluster', 'Disease', 'Seeds'])
    for cluster in filtered_dico_cluster:
        diseases_codes = filtered_dico_cluster[cluster]
        diseases_names = []
        for code in diseases_codes:
            diseases_names.append(dico_code_disease[str(code)])
        for disease, code in zip(diseases_names, diseases_codes):
            seeds = dico_disease_seeds[str(code)]
            new_row = [f"{cluster[8:]}", disease, seeds]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
    df.to_csv(path + f"output_tables/diseases_in_clusters.csv", sep=",", index=False)
    source_file = f"output_tables/diseases_in_clusters.csv"
    destination_folder = path + "../_00_data/"
    shutil.copy(source_file, destination_folder)

analyze_clusters(dico_disease_seeds=dico_disease_seeds, filtered_dico_cluster=filtered_dico_cluster, dico_code_disease=dico_code_disease)