#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 2023

@author: Cécile

Functions to analyze a network interaction file
"""

import networkx as nx
import argparse
import os

# Argparse
parser = argparse.ArgumentParser(
    prog="cluster_communities.py", 
    description="Function to analyze networks"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")

def analyse_network(network: str):
    """Function which analyzes a network
    interactome file by computing some
    network metrics :
       - number of nodes
       - number of edges
       - number of self-loops
       - density

    Args:
        network (str): the name of the interactome file
    """
    # We use the networkx Graph class : undirected graph,
    # wihtout parallel edges
    ntw = nx.read_edgelist(network, create_using = nx.Graph)
    print(f"number of nodes : {nx.number_of_nodes(ntw)}")
    print(f"number of edges :  {nx.number_of_edges(ntw)}")
    print(f"number of self loops : {nx.number_of_selfloops(ntw)}")
    print(f"self loops : {list(nx.nodes_with_selfloops(ntw))}")
    print(f"density : {nx.density(ntw)}")
    """all_degress = ntw.degree()
    average_degree = (sum(all_degress.values()))/(nx.number_of_nodes(ntw))
    print(f"average degree : {average_degree}")"""


print("PPI")
analyse_network(comm_path + "multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv")
print(" ")
print("Pathway")
analyse_network(comm_path + "multiplex/1/reactome_pathways_gene_names_190123.tsv")
print(" ")
print("Complexes")
analyse_network(comm_path + "multiplex/1/Complexes_gene_names_190123.tsv")
print(" ")
print("CoExp")
analyse_network(comm_path + "multiplex/1/Coexpression_310323.tsv")

