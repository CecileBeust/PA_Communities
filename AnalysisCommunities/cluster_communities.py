# Import modules
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import argparse
from analyse_final_communities import get_list_orpha_names, create_dico_disease_seeds

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


# Define the dico of diseases and their seeds + list of ORPHANET identifiers
(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "orpha_codes_PA.txt")

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
        community = comm_path + f"results_{size}_{id}/seeds_{id}.txt"
        # check is path exists = if the RWR has gave an output for this disease
        if os.path.exists(community):
            list_communities.append(community)
        # if not we store the id in a list of ids not analyzed
        else:
            not_analyzed.append(id)
    return (list_communities, not_analyzed)

# Build the list of communities
(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)
print(f"Communities : {communities_10}")
print(" ")
print(f"Diseases not analyzed : {not_analyzed}")
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
print(f"There are {len(list_ids_analyzed)} diseases analyzed")
disease_name_list = get_list_orpha_names(path, "pa_orphanet_diseases.tsv", list_ids_analyzed)
print(" ")
print(f"Diseases names list : {disease_name_list}")

df_disease = pd.read_csv(path + "/pa_orphanet_diseases.tsv", sep="\t", header=None)
disease_names = []
for id in list_ids_analyzed:
    for row in df_disease.iterrows():
        if f"ORPHA:{id}" == row[df_disease.columns[1]][0]:
            disease_names.append(str(row[df_disease.columns[1]][1]))

# Define a list of shorten disease names for better further visualization
shorten_disease_names = ['FG syndrome type 1', 
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
                         'Arterial tortuosity syndrome', 
                         'Trichothiodystrophy', 
                         'Hoyeraal-Hreidarsson syndrome', 
                         'SHORT syndrome', 
                         'Nicolaides-Baraitser syndrome', 
                         'Progeroid syndrome, Petty type', 
                         'Classical Ehlers-Danlos syndrome', 
                         'Vascular Ehlers-Danlos syndrome', 
                         'Wrinkly skin syndrome', 
                         'Moyamoya disease', 
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
                         'Neurologic Waardenburg-Shah syndrome', 
                         'Cataract-intellectual disability-hypogonadism syndrome', 
                         'SOLAMEN Syndrome', 
                         'Cardiofaciocutaneous syndrome', 
                         'Branchioskeletogenital syndrome', 
                         'Branchio-oculo-facial syndrome', 
                         'Transaldolase deficiency', 
                         'Ataxia-telangiectasia'
                         ]

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

# Build similarity matrix of communities based on Jaccard index
(similarity_matrix_10, dico_jaccard) = build_similarity_matrix(communities_10)
print(" ")
print("Similarity matrix based on Jaccard index")
print(similarity_matrix_10)

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
distance_matrix_10 = build_distance_matrix(similarity_matrix_10)

def plot_dendrogram(model, **kwargs):
    # Children of hierarchical clustering
    children = model.children_

    # Distances between each pair of children
    # Since we don't have this information, we can use a uniform one for plotting
    distance = np.arange(children.shape[0])

    # The number of observations contained in each cluster level
    no_of_observations = np.arange(2, children.shape[0]+2)

    # Create linkage matrix and then plot the dendrogram
    linkage_matrix = np.column_stack([children, distance, no_of_observations]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def cluster_dendrogram(matrix: np.ndarray, cutoff: float, size: int):
    """
    Function to plot a hierarchichal clustering dendrogram from a distance matrix

    Args:
        matrix (np.ndarray) : a distance matrix
        cutoff (float) : a value to cut the tree
        size (int) : the number of iterations chose
        to build the communities
    Return:
        None
    """
    # transform distance matrix to flat array, only keep the upper part of the matrix
    array = matrix[np.triu_indices(len(matrix), k=1)]
    # linkage matrix 
    Z = linkage(array, 'average')
    
    # plot dendrogram
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.title("Hierarchichal clustering of disease communities")
    #dn = hierarchy.dendrogram(Z, labels=list_ids_analyzed, leaf_font_size=7, ax=ax)
    #dn = hierarchy.dendrogram(Z, labels=disease_names, leaf_font_size=7, ax=ax, color_threshold=cutoff, orientation='right')
    #dn = dendrogram(Z, labels=shorten_disease_names, leaf_font_size=7, ax=ax, color_threshold=None, orientation='left')
    dn = dendrogram(Z, labels=["HGPS", "WS"], leaf_font_size=7, ax=ax, color_threshold=None, orientation='left')
    # define clusters assignments
    assignments = fcluster(Z, t=cutoff, criterion='distance')
    print(f"Cluster assignements : {assignments}") 
    print(f"Number of clusters : {max(assignments)}")
    print(len(list_ids_analyzed))
    
    # export cluster assignment to dataframe
    cluster_output = pd.DataFrame({'disease': list_ids_analyzed, 'cluster':assignments})
    #cluster_output = pd.DataFrame({'disease': disease_names, 'cluster':assignments})
    # export cluster output to file
    cluster_output.to_csv(f"cluster_output_{size}_{cutoff}.tsv", sep="\t")
    
    # add cutoff line to dendrogram
    plt.axvline(x=cutoff, c='r')
    plt.tight_layout()
    
    # export distance matrix to dataframe for clustermap
    #df = pd.DataFrame(matrix, columns=list_ids_analyzed, index=list_ids_analyzed)
    #df = pd.DataFrame(matrix, columns=disease_names, index=disease_names)
    #df = pd.DataFrame(matrix, columns=shorten_disease_names, index=shorten_disease_names)
    df = pd.DataFrame(matrix, columns=["HGPS", "WS"], index=["HGPS", "WS"])
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
    plt.show()

cluster_dendrogram(distance_matrix_10, 0.7, 10)