from create_config_seed_files import build_seeds_file
import os

dico_seeds = build_seeds_file("orpha_codes_PA.txt")
path = '/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2'

def test_no_seeds(dico_seeds: dict, path: str, list_sizes: list):
    list_id_not_analyzed = []
    for disease in dico_seeds.keys():
        if dico_seeds[disease] == []:
            list_id_not_analyzed.append(disease)
    for size in list_sizes:
        for disease_not_analyzed in list_id_not_analyzed:
            # check this disease has not been analyzed = check result path does not exist
            assert not os.path.exists(f"{path}/results_{size}_{disease_not_analyzed}/seeds_{disease_not_analyzed}.txt")
    return list_id_not_analyzed

#list_id_not_analyzed = test_no_seeds(dico_seeds, path, [100])
list_id_not_analyzed = test_no_seeds(dico_seeds, path, [100])

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
        community_disease_file = f"{path}/results_{size}_{disease}/seeds_{disease}.txt"
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
test_communities(dico_seeds, list_id_not_analyzed, path, 50)

def test_config_files(dico_seeds: dict, not_analyzed: list, path, size):
    disease_analyzed = [disease for disease in dico_seeds.keys() if disease not in not_analyzed]
    for disease in disease_analyzed:
        config_disease = f"{path}/config_{disease}.yml"
        config_community = f"{path}/results_{size}_{disease}/seeds_{disease}.txt"
        with open(config_disease, 'r') as file:
            contents = file.read()
            comp = "seed: seeds_{disease}.txt \
                    self_loops: 0 \
                    r: 0.7 \
                    eta: [1.0] \
                    multiplex: \
                        1: \
                            layers: \
                                - multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv \
                                - multiplex/1/reactome_pathways_gene_names_190123.tsv \
                                - multiplex/1/Coexpression_310323.tsv \
                                - multiplex/1/Complexes_gene_names_190123.tsv \
                                - multiplex/1/Gene_involved_in_diseases_layer_020223.tsv \
                            delta: 0.5 \
                            graph_type: [00, 00, 00, 00, 00] \
                            tau: [0.2, 0.2, 0.2, 0.2, 0.2]"
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
                                - multiplex/1/PPI_HiUnion_LitBM_APID_gene_names_190123.tsv \
                                - multiplex/1/reactome_pathways_gene_names_190123.tsv \
                                - multiplex/1/Coexpression_310323.tsv \
                                - multiplex/1/Complexes_gene_names_190123.tsv \
                                - multiplex/1/Gene_involved_in_diseases_layer_020223.tsv \
                            delta: 0.5 \
                            graph_type: [00, 00, 00, 00, 00] \
                            tau: [0.2, 0.2, 0.2, 0.2, 0.2]"
            assert contents == comp

#test_config_files(dico_seeds, list_id_not_analyzed, path, 100)