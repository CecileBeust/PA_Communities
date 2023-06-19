import os
import pandas as pd
import numpy as np
from itRWR import create_config_seed_files

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

data_folder = os.path.join(os.path.dirname(__file__), '..', 'data')
data_file_path = os.path.join(data_folder, 'orpha_codes_PA.txt')

dico_seeds = create_config_seed_files.build_seeds_file(data_file_path)

def test_no_seeds_diseases(path: str, size: int) -> list:
    """Function to test if the diseases which are not associated to any
    seed are not analyzed during the identification of communities

    Args:
        path (str): path to current Tests folder
        sizes (int): number of iterations used for itRWR
    
    Returns:
        list: a list of diseases not analyzed because not assiciated to
        any seed
    """
    df = pd.read_csv(data_file_path)
    list_id_not_analyzed = list()
    for index, row in df.iterrows():
        if row[2] is np.nan or not pd.isnull(row[2]):
            list_id_not_analyzed.append(row[1])
    for disease_not_analyzed in list_id_not_analyzed:
        # check this disease has not been analyzed = check result path does not exist
        assert not os.path.exists(f"{path}/../IDCommunity/results_{size}_{disease_not_analyzed}/seeds_{disease_not_analyzed}.txt")
    return list_id_not_analyzed

#list_id_not_analyzed = test_no_seeds_diseases(path, 10)

def test_communities(dico_seeds: dict, not_analyzed: list, path, size: int):
    disease_analyzed = [disease for disease in dico_seeds.keys() if disease not in not_analyzed]
    for disease in disease_analyzed:
        seeds_disease_file = f"{path}/seeds_{disease}.txt"
        nb_seeds = 0
        seeds_disease = []
        with open(seeds_disease_file, 'r') as file:
            for line in file:
                nb_seeds += 1
                seeds_disease.append(line.rstrip())
        community_disease_file = f"{path}/../IDCommunity/results_{size}_{disease}/seeds_{disease}.txt"
        community_size = 0
        community_disease_nodes = []
        with open(community_disease_file, 'r') as file:
            for line in file:
                community_size += 1
                community_disease_nodes.append(line.rstrip())
        # check commnunity size
        assert community_size == size + nb_seeds
        # check if seeds are in the community
        for seed in seeds_disease:
            assert seed in community_disease_nodes
        
        
#test_communities(dico_seeds, list_id_not_analyzed, path, 100)

def test_config_files(dico_seeds: dict, not_analyzed: list, path, size):
    disease_analyzed = [disease for disease in dico_seeds.keys() if disease not in not_analyzed]
    for disease in disease_analyzed:
        config_disease = f"{path}/../IDCommunity/config_{disease}.yml"
        config_community = f"{path}../IDCommunity/results_{size}_{disease}/seeds_{disease}.txt"
        with open(config_disease, 'r') as file:
            contents = file.read()
            comp = "seed: seeds_{disease}.txt \
                    self_loops: 0 \
                    r: 0.7 \
                    eta: [1.0] \
                    multiplex: \
                        1: \
                            layers: \
                                - multiplex/1/Coexpression.tsv \
                                - multiplex/1/Complexes.tsv \
                                - multiplex/1/PPI.tsv \
                                - multiplex/1/Pathways.tsv \
                            delta: 0.5 \
                            graph_type: [00, 00, 00, 00] \
                            tau: [0.25, 0.25, 0.25, 0.25]"
            assert contents == comp
        with open(config_community, 'r') as file:
            contents = file.read()
            comp = "seed: results_{size}_{disease}/seeds_{disease}.txt \
                    self_loops: 0 \
                    r: 0.7 \
                    eta: [1.0] \
                    multiplex: \
                        1: \
                            layers: \
                                - multiplex/1/Coexpression.tsv \
                                - multiplex/1/Complexes.tsv \
                                - multiplex/1/PPI.tsv \
                                - multiplex/1/Pathways.tsv \
                            delta: 0.5 \
                            graph_type: [00, 00, 00, 00] \
                            tau: [0.25, 0.25, 0.25, 0.25]"
            assert contents == comp

#test_config_files(dico_seeds, list_id_not_analyzed, path, 10)