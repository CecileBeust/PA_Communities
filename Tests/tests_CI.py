import os
import pandas as pd
import numpy as np
import argparse
from itRWR import create_config_seed_files

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# Argparse
parser = argparse.ArgumentParser(
    prog="cluster_communities.py", 
    description="Test functions"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")

data_folder = os.path.join(os.path.dirname(__file__), '..', '_00_data')
data_file_path = os.path.join(data_folder, 'orpha_codes_PA.txt')

dico_seeds = create_config_seed_files.build_seeds_file(data_file_path)
print(f"Dico seeds: {dico_seeds}")

def test_no_seeds_diseases(size: int) -> list:
    """Function to test if the diseases which are not associated to any
    seed are not analyzed during the identification of communities

    Args:
        sizes (int): number of iterations used for itRWR
    
    Returns:
        list: a list of diseases not analyzed because not assiciated to
        any seed
    """
    df = pd.read_csv(data_file_path, sep="\t", header=0)
    list_id_not_analyzed = list()
    for index, row in df.iterrows():
        if row[1] is np.nan or pd.isnull(row[1]):
            list_id_not_analyzed.append(row[0])
    for disease_not_analyzed in list_id_not_analyzed:
        # check this disease has not been analyzed = check result path does not exist
        assert not os.path.exists(comm_path + f"results_{size}_{disease_not_analyzed}/seeds_{disease_not_analyzed}.txt")
    return list_id_not_analyzed

list_id_not_analyzed = test_no_seeds_diseases(100)

def test_communities(dico_seeds: dict, not_analyzed: list, comm_path: str, size: int) -> None:
    """Function to test the communities: check their size, check if they contain the
    seed(s) node(s)

    Args:
        dico_seeds (dict): dictionary containing diseases ORPHANET codes as keys and
        their causative genes (seeds) as values
        not_analyzed (list): list of ORPHANET codes of diseases not analyzed (because no
        associated gene)
        comm_path (str): path where communities folders are stored (path to _01_Community_Identification)
        size (int): number of iterations used for itRWR, reflecting the sizes of communities
    """
    disease_analyzed = [disease for disease in dico_seeds.keys() if disease not in not_analyzed]
    for disease in disease_analyzed:
        seeds_disease_file = f"{comm_path}/seeds_{disease}.txt"
        nb_seeds = 0
        seeds_disease = []
        # count the number of seeds per disease
        with open(seeds_disease_file, 'r') as file:
            for line in file:
                nb_seeds += 1
                seeds_disease.append(line.rstrip())
        community_disease_file = comm_path + f"results_{size}_{disease}/seeds_{disease}.txt"
        community_size = 0
        community_disease_nodes = []
        # count number of genes in each community
        with open(community_disease_file, 'r') as file:
            for line in file:
                community_size += 1
                community_disease_nodes.append(line.rstrip())
        # check commnunity size
        assert community_size == size + nb_seeds
        # check if seeds are in the community
        for seed in seeds_disease:
            assert seed in community_disease_nodes
        
test_communities(dico_seeds, list_id_not_analyzed, comm_path, 100)

"""def test_config_files(dico_seeds: dict, not_analyzed: list, comm_path: str, size: int):
    disease_analyzed = [disease for disease in dico_seeds.keys() if disease not in not_analyzed]
    for disease in disease_analyzed:
        config_disease = comm_path + f"config_{disease}.yml"
        config_community = comm_path + f"results_{size}_{disease}/config.yml"
        with open(config_disease, 'r') as file:
            contents = file.read()
            print(contents)
            comp = f"seed: seeds_{disease}.txt
            self_loops: 0 + \
            r: 0.7 + \
            eta: [1.0] + \
            multiplex: + \
                1: + \
                    layers: + \
                        - multiplex/1/Coexpression_310323.tsv + \
                        - multiplex/1/Complexes_gene_names_190123.tsv + \
                        - multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv + \
                        - multiplex/1/reactome_pathways_gene_names_190123.tsv + \
                    delta: 0.5 + \
                    graph_type: [00, 00, 00, 00] + \
                    tau: [0.25, 0.25, 0.25, 0.25] "
            assert contents == comp
        with open(config_community, 'r') as file:
            contents = file.read()
            comp = f"seed: results_{size}_{disease}/seeds_{disease}.txt \
                    self_loops: 0 \
                    r: 0.7 \
                    eta: [1.0] \
                    multiplex: \
                        1: \
                            layers: \
                                - multiplex/1/Coexpression_310323.tsv \
                                - multiplex/1/Complexes_gene_names_190123.tsv \
                                - multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv \
                                - multiplex/1/reactome_pathways_gene_names_190123.tsv \
                            delta: 0.5 \
                            graph_type: [00, 00, 00, 00] \
                            tau: [0.25, 0.25, 0.25, 0.25]"
            assert contents == comp

test_config_files(dico_seeds, list_id_not_analyzed, comm_path, 100)"""