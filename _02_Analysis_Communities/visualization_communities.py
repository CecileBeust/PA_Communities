"""
Functions to generate files to visualize
gene communities in Cytoscape
"""

# import modules
import pandas as pd
import networkx as nx
import os
import argparse

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# Argparse
parser = argparse.ArgumentParser(
    prog="visualization_communities.py", 
    description="functions to generate files for the visualization of the communities in Cytoscape"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")

def plot_community(id: int, size: int) -> None:
    """Function to generate a tabulated file
    of nodes in a community and edges
    between them

    Args:
        id (int): the ORPHANET identifier of the
        disease-associated community to analyze
        size (int): number of iterations used to
        build the community with the itRWR algorithm
    """
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    community = comm_path + f"results_{size}_{id}/seeds_{id}.txt"
    nodes = []
    with open(community, 'r') as community_file:
        for node in community_file:
            nodes.append(node.rstrip())
    df = pd.DataFrame(columns=['node1', 'node2', 'provenance'])
    # extract subnetwork with commmunity nodes for each layer
    nodes_ppi = ppi.subgraph(nodes)
    edges_ppi = nodes_ppi.edges()
    nodes_pathways = pathways.subgraph(nodes)
    edges_pathways = nodes_pathways.edges()
    nodes_coexp = coexp.subgraph(nodes)
    edges_coexp = nodes_coexp.edges()
    nodes_complexes = complexes.subgraph(nodes)
    edges_complexes = nodes_complexes.edges()
    int_ppi = []
    int_pathways = []
    int_coexp = []
    int_complexes = []
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
    df.to_csv(f"output_visualization/{id}_community_{size}.tsv", sep="\t")
    print(df)

# Here we want to generate the files for the visualisation of the communities
# of Hutchinson-Gilford Progeria Syndrome (ORPHANET code 740) 

# Change the ORPHANET code depending on the community to analyze
plot_community(740, 100)

def plot_community_with_neighbors(id: int, size: int) -> None:
    """Function to generate a tabulated file of the genes in 
    a community and their neighbors, with the edeges between
    them

    Args:
        id (int): _description_
        size (int): _description_
    """
    ppi = nx.read_edgelist(comm_path + "multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + "multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + "multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + "multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
    community = comm_path + f"results_{size}_{id}/seeds_{id}.txt"
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
    print(counter)
    print(df)
    df.to_csv(f"output_visualization/{id}_community_{size}_with_neighbors.tsv", sep="\t")


def plot_several_communities(list_id: list, size: int) -> set:
    """Function to generate a tsv file for the visualization
    of several communities together in Cytoscape

    Args:
        list_id (list): list of ORPHANET codes of diseases
        to analyze
        size (int): number of iterations used for itRWR, reflecting
        the size of communities

    Returns:
        set: set of nodes in the communities analyzed
    """
    ppi = nx.read_edgelist(comm_path + 'multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv', create_using = nx.Graph)
    pathways = nx.read_edgelist(comm_path + 'multiplex/1/reactome_pathways_gene_names_190123.tsv', create_using = nx.Graph)
    coexp = nx.read_edgelist(comm_path + 'multiplex/1/Coexpression_310323.tsv', create_using = nx.Graph)
    complexes = nx.read_edgelist(comm_path + 'multiplex/1/Complexes_gene_names_190123.tsv', create_using = nx.Graph)
    df = pd.DataFrame(columns=['node1', 'node2', 'provenance'])
    #extract all nodes in the 3 communities and the associated subnetwork
    nodes_all_comm = set()
    for id in list_id :
        community = comm_path + f'/results_{size}_{id}/seeds_{id}.txt'
        nodes = []
        with open(community, 'r') as community_file:
            for node in community_file:
                nodes.append(node.rstrip())
                nodes_all_comm.add(node.rstrip())
    nodes_ppi_all = ppi.subgraph(list(nodes_all_comm))
    edges_ppi_all = nodes_ppi_all.edges()
    nodes_pathways_all = pathways.subgraph(list(nodes_all_comm))
    edges_pathways_all = nodes_pathways_all.edges()
    nodes_coexp_all = coexp.subgraph(list(nodes_all_comm))
    edges_coexp_all = nodes_coexp_all.edges()
    nodes_complexes_all = complexes.subgraph(list(nodes_all_comm))
    edges_complexes_all = nodes_complexes_all.edges()
    int_ppi_all = []
    int_pathways_all = []
    int_coexp_all = []
    int_complexes_all = []
    counter = 0
    for interaction in edges_ppi_all:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_ppi_all and not reverse in int_ppi_all:
            new_row = [interaction[0], interaction[1], 'ppi']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_ppi_all.append(forward)
    for interaction in edges_pathways_all:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_pathways_all and not reverse in int_pathways_all:
            new_row = [interaction[0], interaction[1], 'pathways']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_pathways_all.append(forward)
    for interaction in edges_coexp_all:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_coexp_all and not reverse in int_coexp_all:
            new_row = [interaction[0], interaction[1], 'coexp']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_coexp_all.append(forward)
    for interaction in edges_complexes_all:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_complexes_all and not reverse in int_complexes_all:
            new_row = [interaction[0], interaction[1], 'complexes']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_complexes_all.append(forward)
    df.to_csv(path + f'output_visualization/3communities_{size}.tsv', sep="\t", index=False)
    return nodes_all_comm 
    

# Use this command to generate the tsv file for the representation of the three communities of
# Hutchinson-Gilford Progeria Syndrome (ORPHANET code 740), Ataxia telangiectasia (ORPHANET code 100)
# and Classical Ehlers-Danlos syndrome (OPRHANET code 287)
nodes_all_comm = plot_several_communities(['100', '740', '287'], 100)