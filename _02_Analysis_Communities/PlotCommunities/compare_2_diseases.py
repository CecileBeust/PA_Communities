import pandas as pd
import networkx as nx
import numpy as np

def plot_community(id, size):
    ppi = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist("/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)
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
    df.to_csv(f"{id}_community_{size}.tsv", sep="\t", index=None)
    print(counter)
    print(df)

# werner syndrome
#plot_community(902, 100)

# classical ehlers danlos syndrome
#plot_community(287, 100)

def create_aggregate_ntw():
    ppi = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv"
    pathways = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/reactome_pathways_gene_names_190123.tsv"
    coexp = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Coexpression_310323.tsv"
    complexes = "/home/cbeust/Landscape_PA/CommunityIdentification/RESULTS_CI/multiplex/1/Complexes_gene_names_190123.tsv"
    diseases = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv"
    edgelists = [ppi, pathways, coexp, complexes, diseases]
    agg_graph = nx.Graph()
    for edgelist in edgelists:
        G = nx.read_edgelist(edgelist, create_using = nx.Graph)
        agg_graph.add_edges_from(G.edges())
    # Remove self-loops from the aggregated graph
    self_loops = nx.selfloop_edges(agg_graph, data=True) 
    agg_graph.remove_edges_from(self_loops)
    # Remove duplicate edges from the aggregated graph
    # agg_graph = nx.unique_edges(agg_graph)
    print(f"The aggregated network has {agg_graph.number_of_nodes()} nodes and {agg_graph.number_of_edges()} edges")
    nx.write_edgelist(agg_graph, "aggregated_network.tsv", delimiter="\t")

#create_aggregate_ntw()

"""with open("aggregated_network.tsv", 'r') as file:
    edges = list()
    for line in file:
        gene1 = line.split("\t")[0].rsplit()
        gene2 = line.split("\t")[1].rsplit()
        interaction = (gene1, gene2)
        reverse = (gene2, gene1)
        if interaction in edges or reverse in edges:
            print(interaction)
        else:
            edges.append(interaction)
print(len(interaction))"""

def create_node_list_comm(id1, id2, size):
    community1 = f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/results_{size}_{id1}/seeds_{id1}.txt"
    community2 = f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/results_{size}_{id2}/seeds_{id2}.txt"
    nodes = set()
    with open(community1, 'r') as file:
        for line in file:
            nodes.add(line.rstrip())
    with open(community2, 'r') as file:
        for line in file:
            nodes.add(line.rstrip())
    with open(f"{id1}_{id2}_nodes.txt", 'w') as out:
        for node in list(nodes):
            out.write(str(node) + "\n")

#create_node_list_comm(902, 287, 100)

def add_node_property(agg_ntw, id1, id2, size):
    df = pd.read_csv(agg_ntw, sep="\t", names=['node1', 'node2', 'other'])
    df['provenance'] = 0
    print(df)
    comm1 = pd.read_csv(f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/Analysis_Communities_V2/{id1}_community_{size}.tsv", sep="\t")
    comm1.drop(3)
    comm1_int = set()
    comm2 = pd.read_csv(f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/Analysis_Communities_V2/{id2}_community_{size}.tsv", sep="\t")
    comm2.drop(3)
    comm2_int = set()
    for index, row in comm1.iterrows():
        interaction = row["node1"] + "\t" + row["node2"] + "\n"
        comm1_int.add(interaction)
    for index, row in comm2.iterrows():
        interaction = row["node1"] + "\t" + row["node2"] + "\n"
        comm2_int.add(interaction)
    both = comm1_int.intersection(comm2_int)
    print(comm1_int)
    i = 0
    for index, row in df.iterrows():
        interaction = row[0] + "\t" + row[1] + "\n"
        reverse = row[1] + "\t" + row[0] + "\n"
        if interaction in list(comm1_int) or reverse in list(comm1_int):
            if interaction not in list(both):
                df._set_value(i, "provenance", id1)
        if interaction in list(comm2_int) or reverse in list(comm2_int):
            if interaction not in list(both):
                df._set_value(i, "provenance", id2)
        if interaction in list(both) or reverse in list(both):
            df._set_value(i, "provenance", "both")
        i += 1
    print(df)
    df.to_csv("aggregated_prop_comm.tsv", sep="\t", index=None)

add_node_property("aggregated_network.tsv", 902, 287, 100)