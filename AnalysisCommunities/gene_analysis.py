from analyse_final_communities import create_dico_disease_seeds, build_communities_list, build_similarity_matrix, jaccard_index
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import collections
import os

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

comm_path = "/home/cbeust/Landscape_PA/Github_Codes/PA_Communities/IDCommunity"

(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "orpha_codes_PA.txt")
(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]
print(list_ids_analyzed)
similarity_matrix_10 = build_similarity_matrix(communities_10)[0]
print(similarity_matrix_10)

df_disease = pd.read_csv("pa_orphanet_diseases.tsv", sep="\t", header=None)
disease_names = []
for id in list_ids_analyzed:
    for row in df_disease.iterrows():
        if f"ORPHA:{id}" == row[df_disease.columns[1]][0]:
            disease_names.append(str(row[df_disease.columns[1]][1]))

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


def plot_clustermap(sim_matrix: np.ndarray, value) -> None:
    """Function to build a clustermap from a 
    similarity matrix

    Args:
        sim_matrix (np.ndarray): a similarity matrix
        value (_type_): number of iterations to build the 
        communities --> for the title of the clustermap
    Return:
        None
    """
    df = pd.DataFrame(sim_matrix, columns=list_ids_analyzed, index=list_ids_analyzed)
    print(df)
    clust = sns.clustermap(
        df,
        annot=False,
        metric="euclidean",
        method="average",
        cmap="magma_r",
        row_cluster=True,
        col_cluster=True,
        figsize=(10,10),
    )
    # set titles and legend
    clust.ax_heatmap.set_ylabel("Premature Aging Diseases Orphanet Codes")
    if value == 30:
        clust.fig.suptitle("Jaccard index clustermap of size 30 + seeds communities identified by RWR using Premature Aging disease genes as seeds")
    elif value == 50:
        clust.fig.suptitle("Jaccard index clustermap of size 50 + seeds communities identified by RWR using Premature Aging disease genes as seeds")
    elif value == 100:
        clust.fig.suptitle("Jaccard index clustermap of size 100 + seeds communities identified by RWR using Premature Aging disease genes as seeds")
    elif value == 0:
        pass

    plt.show()

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

def analyse_genes_in_modules(communities, list_ids_analyzed, dico_seeds):
    dico_gene_comm = {}
    genes_total = []
    genes_wo_seeds = []
    genes_in_several_modules = []
    seeds = set()
    for values in dico_seeds.values():
        for seed in values:
            seeds.add(seed)
    for module, id in zip(communities, list_ids_analyzed):
        with open(module, 'r') as file:
            for line in file:
                gene = line.rstrip()
                # if the gene has not been registered in the lits of total genes
                if not gene in genes_total:
                    # add it to the list of genes and register it in the dico
                    genes_total.append(gene)
                    dico_gene_comm[gene] = [module]
                # if the gene has been registered
                else:
                    dico_gene_comm[gene] += [module]
                    # if the gene has not been registered in the list of genes appearing in several modules
                    if gene not in genes_in_several_modules:
                        # add it to the list of genes involved in several modules
                        genes_in_several_modules.append(gene)
                # if the gene is not a seed
                if not gene in seeds:
                    if not gene in genes_wo_seeds:
                        genes_wo_seeds.append(gene)
    #print(dico_gene_comm)
    return genes_total, genes_wo_seeds, dico_gene_comm

"""
analyse_genes_in_modules(communities_30, list_ids_analyzed, dico_disease_seeds)
analyse_genes_in_modules(communities_50, list_ids_analyzed, dico_disease_seeds)
"""
"""genes_100 = analyse_genes_in_modules(communities_100, list_ids_analyzed, dico_disease_seeds)[0]
print(len(genes_100))
genes_100_wo_seeds = analyse_genes_in_modules(communities_100, list_ids_analyzed, dico_disease_seeds)[1]
print(len(genes_100_wo_seeds))
dico_genes_comm_100 = analyse_genes_in_modules(communities_100, list_ids_analyzed, dico_disease_seeds)[2]
print(len(dico_genes_comm_100))

df_100 = pd.DataFrame(np.zeros((1770, 2)))
df_100.columns = ['Gene name', 'Number of PA communities memberships']
i = 0
for gene in dico_genes_comm_100.keys():
    if gene in genes_100_wo_seeds:
        df_100._set_value(i, 'Gene name', gene)
        nb_comm = len(dico_genes_comm_100[gene])
        df_100._set_value(i, 'Number of PA communities memberships', nb_comm)
        i += 1
#df_100.sort_values(by=['Number of PA communities memberships'], ascending=False)
print(df_100)
df_100.to_csv("genes_comm_100.tsv", sep="\t", index=None)"""

def generate_excel_genes(dico_gene_comm):
    df = pd.DataFrame(np.zeros((1903, 9)))
    df.columns = ['Gene symbol', 'Nb of communities', 'degree PPI ntw', 'degree Pathways ntw', 'degree Coexp ntw', 'degree Complexes ntw', 'degree Disease ntw', 'max degree', 'sum degrees']
    i = 0
    ppi = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    diseases = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv", create_using = nx.Graph)
    for gene in dico_gene_comm:
        df._set_value(i, 'Gene symbol', gene)
        nb_comm = 0
        comm_names = []
        for comm in dico_gene_comm[gene]:
            nb_comm += 1
            comm_names.append(comm)
            df._set_value(i, 'Nb of communities', nb_comm)
            #df._set_value(i, 'Communities', comm_names)
        all_degs = []
        for elt in ppi.degree([gene]):
            deg_ppi = elt[1]
            all_degs.append(deg_ppi)
            df._set_value(i, 'degree PPI ntw', deg_ppi)
        for elt in pathways.degree([gene]):
            deg_pathways = elt[1]
            all_degs.append(deg_pathways)
            df._set_value(i, 'degree Pathways ntw', deg_pathways)
        for elt in coexp.degree([gene]):
            deg_coexp = elt[1]
            all_degs.append(deg_coexp)
            df._set_value(i, 'degree Coexp ntw', deg_coexp)
        for elt in complexes.degree([gene]):
            deg_complexes = elt[1]
            all_degs.append(deg_complexes)
            df._set_value(i, 'degree Complexes ntw', deg_complexes)
        for elt in diseases.degree([gene]):
            deg_diseases = elt[1]
            all_degs.append(deg_diseases)
            df._set_value(i, 'degree Disease ntw', deg_diseases)
        df._set_value(i, 'max degree', max(all_degs))
        df._set_value(i, 'sum degrees', sum(all_degs))
        i += 1
    print(df)
    # Pearson correlation test
    comm = df['Nb of communities']
    max_deg = df['max degree']
    sum_deg = df['sum degrees']
    corr_max, p_value_max = pearsonr(comm, max_deg)
    print('Correlation Nb comm / max deg:', corr_max)
    print('P-value:', p_value_max)
    corr_sum, p_value_sum = pearsonr(comm, sum_deg)
    print('Correlation Nb comm / sum deg:', corr_sum)
    print('P-value:', p_value_sum)
    df.to_csv("gene_in_communities.tsv", sep="\t")

#generate_excel_genes(dico_genes_comm)

def plot_community_with_neighbors(id, size):
    ppi = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/Coexpresssion_th05_020223_final.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    diseases = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv", create_using = nx.Graph)
    community = f"/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/results_{size}_{id}/seeds_{id}.txt"
    nodes = []
    with open(community, 'r') as community_file:
        for node in community_file:
            nodes.append(node.rstrip())
    df = pd.DataFrame(columns=['node1', 'node2', 'provenance'])
    counter = 0
    int_ppi = []
    int_pathways = []
    int_coexp = []
    int_complexes = []
    int_disease = []
    for node in nodes:
        for interaction in ppi.edges(node):
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_ppi and not reverse in int_ppi:
                new_row = [interaction[0], interaction[1], 'ppi']
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_ppi.append(forward)
        for interaction in pathways.edges(node):
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_pathways and not reverse in int_pathways:
                new_row = [interaction[0], interaction[1], 'pathways']
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_pathways.append(forward)
        for interaction in coexp.edges(node):
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_coexp and not reverse in int_coexp:
                new_row = [interaction[0], interaction[1], 'coexp']
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_coexp.append(forward)
        for interaction in complexes.edges(node):
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_complexes and not reverse in int_complexes:
                new_row = [interaction[0], interaction[1], 'complexes']
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_complexes.append(forward)
        for interaction in diseases.edges(node):
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_disease and not reverse in int_disease:
                new_row = [interaction[0], interaction[1], 'diseases']
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_disease.append(forward)
    print(counter)
    print(df)
    df.to_csv(f"{id}_community_{size}.tsv", sep="\t")

def plot_community(id, size):
    ppi = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    diseases = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv", create_using = nx.Graph)
    community = f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/results_{size}_{id}/seeds_{id}.txt"
    nodes = []
    with open(community, 'r') as community_file:
        for node in community_file:
            nodes.append(node.rstrip())
    df = pd.DataFrame(columns=['node1', 'node2', 'provenance'])
    # extract subnetwork with commmunity nodes for each layer
    nodes_ppi = ppi.subgraph(nodes)
    edges_ppi = nodes_ppi.edges()
    print(edges_ppi)
    nodes_pathways = pathways.subgraph(nodes)
    edges_pathways = nodes_pathways.edges()
    print(" ")
    print(edges_pathways)
    nodes_coexp = coexp.subgraph(nodes)
    edges_coexp = nodes_coexp.edges()
    print(" ")
    print(edges_coexp)
    nodes_complexes = complexes.subgraph(nodes)
    edges_complexes = nodes_complexes.edges()
    print(" ")
    print(edges_complexes)
    nodes_diseases = diseases.subgraph(nodes)
    edges_diseases = nodes_diseases.edges()
    print(" ")
    print(edges_diseases)
    # edit dataframe
    int_ppi = []
    int_pathways = []
    int_coexp = []
    int_complexes = []
    int_disease = []
    counter = 0
    for interaction in edges_ppi:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_ppi and not reverse in int_ppi:
            new_row = [interaction[0], interaction[1], 'ppi']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_ppi.append(forward)
    for interaction in edges_pathways:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_pathways and not reverse in int_pathways:
            new_row = [interaction[0], interaction[1], 'pathways']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_pathways.append(forward)
    for interaction in edges_coexp:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_coexp and not reverse in int_coexp:
            new_row = [interaction[0], interaction[1], 'coexp']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_coexp.append(forward)
    for interaction in edges_complexes:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_complexes and not reverse in int_complexes:
            new_row = [interaction[0], interaction[1], 'complexes']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_complexes.append(forward)
    for interaction in edges_diseases:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_disease and not reverse in int_disease:
            new_row = [interaction[0], interaction[1], 'diseases']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_disease.append(forward)
    df.to_csv(f"{id}_community_{size}.tsv", sep="\t")
    print(counter)
    print(df)


#plot_community(79474, 50)
