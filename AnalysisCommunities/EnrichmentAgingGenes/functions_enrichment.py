import networkx as nx
import pandas as pd
import numpy as np
import scipy.stats as stats
path = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3"


#############################################
###### LOAD AND PROCESS DATA ################
#############################################

def load_networks():
    # load networks
    ppi = nx.read_edgelist(path + "/multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv", create_using = nx.Graph)
    pathways = nx.read_edgelist(path + "/multiplex/1/reactome_pathways_gene_names_190123.tsv", create_using = nx.Graph)
    coexp = nx.read_edgelist(path + "/multiplex/1/Coexpression_310323.tsv", create_using = nx.Graph)
    complexes = nx.read_edgelist(path + "/multiplex/1/Complexes_gene_names_190123.tsv", create_using = nx.Graph)

    # get noeds in networks
    ppi_nodes = ppi.nodes()
    pat_nodes = pathways.nodes()
    coexp_nodes = coexp.nodes()
    complexes_nodes = complexes.nodes()
    #all_nodes = list(set(list(ppi_nodes) + list(pat_nodes) + list(coexp_nodes) + list(complexes_nodes) + list(diseases_nodes)))
    all_nodes = np.unique(list(ppi_nodes) + list(pat_nodes) + list(coexp_nodes) + list(complexes_nodes))
    print(f"{len(all_nodes)} nodes in the multiplex network")
    return all_nodes

def extract_seeds(file, list_ids_analyzed):
    orpha = pd.read_table(file)
    print(orpha)
    seeds = set()
    for index, row in orpha.iterrows():
        if str(row[0]) in list_ids_analyzed:
            seed = row[1].split(",")
            print(seed)
            for s in seed:
                seeds.add(s)
    return seeds

def load_geneage(file, seeds, nodes_ntw):
    # Load and filter Gene Age database
    df = pd.read_csv(file, sep=",")
    # remove genes without direct association with aging
    filtered_df = df.loc[(df['why'] != "upstream") & (df['why'] != "downstream") & (df['why'] != "putative") & (df['why'] != "functional") & (df['why'] != "downstream,putative") & (df['why'] != "functional,downstream") & (df['why'] != "functional,putative") & (df['why'] != "functional,upstream") & (df['why'] != "upstream,putative")]
    genes_aging = filtered_df['symbol'].to_list()
    print(f"{len(genes_aging)} aging genes")
    # remove seeds from aging genes to avoid overfitting
    genes_aging_wo_seeds = []
    for gene in genes_aging:
        if not gene in seeds:
            genes_aging_wo_seeds.append(gene)
    print(f"{len(genes_aging_wo_seeds)} aging genes after removing seeds genes")
    # check if aging genes are present in networks
    # remove aging gene not present in network
    for gene in genes_aging_wo_seeds:
        if gene not in nodes_ntw:
            genes_aging_wo_seeds.remove(str(gene))
    print(f"{len(genes_aging_wo_seeds)} aging genes without seeds in network")
    return genes_aging_wo_seeds

def extract_genes_from_comm(size: int, list_id_analyzed: list):
    dico_comm_nodes = {}
    for id in list_id_analyzed:
        community = path + f"/results_{size}_{id}/seeds_{id}.txt"
        nodes_comm = []
        with open(community, 'r') as fi:
            for line in fi:
                gene = line.rsplit()
                nodes_comm.append(gene)
                dico_comm_nodes[id] = []
        for node in nodes_comm:
            dico_comm_nodes[id] += node
    return dico_comm_nodes

def extract_genes_from_cluster(size: int, cluster_assignement: str):
    dico_clusters = {}
    filtered_dico = {}
    df = pd.read_csv(cluster_assignement, sep="\t", header=None, index_col=None)
    #df.drop(columns=df.columns[0], axis=1, inplace=True)
    df = df.reset_index()
    print(df)
    for index, row in df.iterrows():
        if not row[2] in dico_clusters.keys():
            dico_clusters[row[2]] = [row[1]]
        else:
            dico_clusters[row[2]] += [row[1]]
    for cluster in dico_clusters:
        if len(dico_clusters[cluster]) >= 3:
            filtered_dico[cluster] = dico_clusters[cluster]
    dico_clusters_nodes = {}
    for cluster in filtered_dico:
        dico_clusters_nodes[cluster] = []
        for comm in filtered_dico[cluster]:
            community = path + f"/results_{size}_{comm}/seeds_{comm}.txt"
            with open(community, 'r') as file:
                for line in file:
                    gene = line.rstrip()
                    if not gene in dico_clusters_nodes[cluster]:
                        dico_clusters_nodes[cluster] += [gene]
    return dico_clusters_nodes

##################################################
######## ENRICHMENT FUNCTIONS ####################
##################################################

# Fischer's exact test : compare distributions of genes_comm and genes_aging_wo_seeds
def fisher(list1, list2, gene_pool):
    common_genes = set(list1).intersection(set(list2))
    print(f"{len(common_genes)} common genes")
    set1_unique = set(list1) - set(list2)
    print(f"{len(set1_unique)} genes in genes_comm but not in genes_aging_we_seeds")
    set2_unique = set(list2) - set(list1)
    print(f"{len(set2_unique)} genes in genes_aging_wo_seeds but not in genes_comm")

    not_set2 = gene_pool - len(set2_unique) - len(common_genes)
    #neither = gene_pool - len(common_genes) - len(set1_unique) - len(set2_unique)
    contingency_table = [
        [len(common_genes), len(set2_unique)],
        [len(set1_unique), not_set2]
    ]

    odds_ratio, p_value = stats.fisher_exact(contingency_table)

    print(f'Odds Ratio: {odds_ratio:.10f}, p-value: {p_value:.10f}')
    return p_value


def hypergeome(list1, list2, gene_pool):
    # Define the size of the gene pool (total number of genes)
    gene_pool_size = gene_pool

    # Define the number of genes in the two lists
    list1_size = len(list1)
    list2_size = len(list2)

    # Determine the number of genes that are common to both lists
    common_genes = set(list1).intersection(list2)
    common_genes_size = len(common_genes)

    # Define the number of genes to randomly sample from the gene pool
    sample_size = list1_size

    # Calculate the p-value using a hypergeometric test
    p_value = stats.hypergeom.sf(common_genes_size-1, gene_pool_size, list2_size, sample_size)

    print(f'p-value: {p_value:.4f}')
    return p_value


def overlap_sets(setA, setB):
    """
    Accepts to lists
    M is the population size (previously N)
    n is the number of successes in the population 
    N is the sample size (previously n)
    x is still the number of drawn “successes”
    """

    M= 18629 #total number of genes in the genome 
    n= len(setA)
    N= len(setB)
    x= len(setA.intersection(setB))
 

    print('p-value <= ' + str(x) + ': ' + str(stats.hypergeom.cdf(x, M, n, N)))
    print('p-value >= ' + str(x)  + ': ' + str(stats.hypergeom.sf(x-1, M, n, N)))
    pass 
