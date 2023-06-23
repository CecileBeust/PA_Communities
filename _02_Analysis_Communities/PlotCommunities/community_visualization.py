import networkx as nx
import pandas as pd

path = '/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3'

def plot_communities(list_id, size):
    ppi = nx.read_edgelist(path + '/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv', create_using = nx.Graph)
    pathways = nx.read_edgelist(path + '/multiplex/1/reactome_pathways_gene_names_190123.tsv', create_using = nx.Graph)
    coexp = nx.read_edgelist(path + '/multiplex/1/Coexpression_310323.tsv', create_using = nx.Graph)
    complexes = nx.read_edgelist(path + '/multiplex/1/Complexes_gene_names_190123.tsv', create_using = nx.Graph)
    #diseases = nx.read_edgelist(path +  '/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv', create_using = nx.Graph)
    df = pd.DataFrame(columns=['node1', 'node2', 'provenance'])
    #extract all nodes in the 3 communities and the associated subnetwork
    nodes_all_comm = set()
    for id in list_id :
        community = path + f'/results_{size}_{id}/seeds_{id}.txt'
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
    """nodes_diseases_all = diseases.subgraph(list(nodes_all_comm))
    edges_diseases_all = nodes_diseases_all.edges()"""
    int_ppi_all = []
    int_pathways_all = []
    int_coexp_all = []
    int_complexes_all = []
    #int_disease_all = []
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
    """for interaction in edges_diseases_all:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_disease_all and not reverse in int_disease_all:
            new_row = [interaction[0], interaction[1], 'diseases']
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_disease_all.append(forward)"""
    """for id in list_id:
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
                new_row = [interaction[0], interaction[1], 'ppi', id, id, id]
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_ppi.append(forward)
        for interaction in edges_pathways:
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_pathways and not reverse in int_pathways:
                new_row = [interaction[0], interaction[1], 'pathways', id, id, id]
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_pathways.append(forward)
        for interaction in edges_coexp:
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_coexp and not reverse in int_coexp:
                new_row = [interaction[0], interaction[1], 'coexp', id, id, id]
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_coexp.append(forward)
        for interaction in edges_complexes:
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_complexes and not reverse in int_complexes:
                new_row = [interaction[0], interaction[1], 'complexes', id, id, id]
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_complexes.append(forward)
        for interaction in edges_diseases:
            forward = (interaction[0], interaction[1])
            reverse = (interaction[1], interaction[0])
            if not forward in int_disease and not reverse in int_disease:
                new_row = [interaction[0], interaction[1], 'diseases', id, id, id]
                df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
                counter += 1
                int_disease.append(forward)"""
    df.to_csv(path + f'/Analysis_Communities_V3/VisualisationComm/3communities_{size}.tsv', sep="\t", index=False)
    return nodes_all_comm 
    

nodes_all_comm = plot_communities(['100', '740', '287'], 100)

def plot_comm(id, size):
    ppi = nx.read_edgelist(path + '/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv', create_using = nx.Graph)
    pathways = nx.read_edgelist(path + '/multiplex/1/reactome_pathways_gene_names_190123.tsv', create_using = nx.Graph)
    coexp = nx.read_edgelist(path + '/multiplex/1/Coexpression_310323.tsv', create_using = nx.Graph)
    complexes = nx.read_edgelist(path + '/multiplex/1/Complexes_gene_names_190123.tsv', create_using = nx.Graph)
    #diseases = nx.read_edgelist(path +  '/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv', create_using = nx.Graph)
    df = pd.DataFrame(columns=['node1', 'node2', 'provenanceedge', "disease"])
    community = path + f'/results_{size}_{id}/seeds_{id}.txt'
    nodes = []
    with open(community, 'r') as community_file:
        for node in community_file:
            nodes.append(node.rstrip())
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
    """nodes_diseases = diseases.subgraph(nodes)
    edges_diseases = nodes_diseases.edges()
    print(" ")
    print(edges_diseases)"""
    # edit dataframe
    int_ppi = []
    int_pathways = []
    int_coexp = []
    int_complexes = []
    #int_disease = []
    counter = 0
    for interaction in edges_ppi:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_ppi and not reverse in int_ppi:
            new_row = [interaction[0], interaction[1], 'ppi', id]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_ppi.append(forward)
    for interaction in edges_pathways:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_pathways and not reverse in int_pathways:
            new_row = [interaction[0], interaction[1], 'pathways', id]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_pathways.append(forward)
    for interaction in edges_coexp:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_coexp and not reverse in int_coexp:
            new_row = [interaction[0], interaction[1], 'coexp', id]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_coexp.append(forward)
    for interaction in edges_complexes:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_complexes and not reverse in int_complexes:
            new_row = [interaction[0], interaction[1], 'complexes', id]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_complexes.append(forward)
    """for interaction in edges_diseases:
        forward = (interaction[0], interaction[1])
        reverse = (interaction[1], interaction[0])
        if not forward in int_disease and not reverse in int_disease:
            new_row = [interaction[0], interaction[1], 'diseases', id]
            df = df.append(dict(zip(df.columns, new_row)), ignore_index=True)
            counter += 1
            int_disease.append(forward)"""
    print(df)
    df.to_csv(path + f'/Analysis_Communities_V3/VisualisationComm/{id}_comm_{size}.tsv', sep="\t", index=False)
    return nodes

nodes_100 = plot_comm("100", 100)
nodes_740 = plot_comm("740", 100)
nodes_287 = plot_comm("287", 100)

node_table = pd.DataFrame(columns=['node', 'provenance'])
nodes_seen = []
for node in nodes_all_comm:
    if node not in nodes_seen:
        if node in nodes_100 and node in nodes_740:
            new_row = [node, "AT;HGPS"]
        elif node in nodes_100 and node in nodes_287:
            new_row = [node, "AT;ED"]
        elif node in nodes_740 and node in nodes_287:
            new_row = [node, "HGPS;ED"]
        elif node in nodes_100:
            new_row = [node, "AT"]
        elif node in nodes_740:
            new_row = [node, "HGPS"]
        elif node in nodes_287:
            new_row = [node, "ED"]
        node_table = node_table.append(dict(zip(node_table.columns, new_row)), ignore_index=True)
        nodes_seen.append(node)
print(node_table)
node_table.to_csv(path + '/Analysis_Communities_V3/VisualisationComm/node_table_3comm.tsv', sep="\t", index=False)


def plot_comm_stats(id, size):
    ppi = nx.read_edgelist(path + '/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv', create_using = nx.Graph)
    pathways = nx.read_edgelist(path + '/multiplex/1/reactome_pathways_gene_names_190123.tsv', create_using = nx.Graph)
    coexp = nx.read_edgelist(path + '/multiplex/1/Coexpression_310323.tsv', create_using = nx.Graph)
    complexes = nx.read_edgelist(path + '/multiplex/1/Complexes_gene_names_190123.tsv', create_using = nx.Graph)
    #diseases = nx.read_edgelist(path +  '/multiplex/1/Gene_involved_in_diseases_layer_020223.tsv', create_using = nx.Graph)
    df = pd.DataFrame(columns=['node1', 'node2', 'provenanceedge', "disease"])
    community = path + f'/results_{size}_{id}/seeds_{id}.txt'
    nodes = []
    with open(community, 'r') as community_file:
        for node in community_file:
            nodes.append(node.rstrip())
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