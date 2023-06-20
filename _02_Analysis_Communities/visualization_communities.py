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
    df.to_csv(f"{id}_community_{size}.tsv", sep="\t")
    print(df)

plot_community(740, 10)
plot_community(902, 10)

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
    df.to_csv(f"{id}_community_{size}_with_neighbors.tsv", sep="\t")

plot_community_with_neighbors(740, 100)
plot_community_with_neighbors(902, 100) 