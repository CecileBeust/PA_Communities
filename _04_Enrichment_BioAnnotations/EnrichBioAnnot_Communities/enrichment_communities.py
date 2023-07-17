import pandas as pd
from gprofiler import GProfiler
import os
import ast
import argparse
import sys
from pathlib import Path

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../..')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, build_communities_list, load_networks

data_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
mapping_hugo_file = os.path.join(data_folder, 'mapping_HGNC_Ensembl.txt')

# Argparse
parser = argparse.ArgumentParser(
    prog="enrichment_communities.py", 
    description="functions to perform enrichment of gene communities"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")


(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]

all_nodes_HUGO = list(load_networks(comm_path=comm_path))

def create_mapping_dict(mapping_file: str):
    mapping_dict = dict()
    df = pd.read_table(mapping_file, header=0)
    for index, row in df.iterrows():
        if str(row[0]) not in mapping_dict.keys():
            mapping_dict[str(row[0])] = str(row[1])
    return mapping_dict

mapping_dict = create_mapping_dict(mapping_file=mapping_hugo_file)

def map_all_nodes_mx(all_nodes: list, mapping_dict: dict):
    all_nodes_ensembl = list()
    for node in all_nodes:
        if node in mapping_dict:
            ensembl_node = mapping_dict[node]
            all_nodes_ensembl.append(ensembl_node)
        else:
            all_nodes_ensembl.append(node)
    return all_nodes_ensembl

all_nodes_ensembl = map_all_nodes_mx(all_nodes=all_nodes_HUGO, mapping_dict=mapping_dict)
print(len(all_nodes_ensembl))

pa_diseases = pd.read_csv(orpha_names, sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

def enrichment_communities(list_comm, list_ids_analyzed, size):
    #os.mkdir('output_tables')
    print(len(list_comm))
    for comm, id in zip(list_comm, list_ids_analyzed):
        with open(comm, 'r') as file:
            genes = []
            for line in file:
                genes += line.rsplit()
            gp = GProfiler(user_agent='https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/gost/profile/', return_dataframe=True)
            enrich = gp.profile(organism='hsapiens', query=genes, sources=["GO:BP", "GO:CC", "REAC"], user_threshold=0.05, domain_scope='custom', significance_threshold_method='fdr', no_evidences=False, background=all_nodes_ensembl)
            print(enrich)
            enrich.to_csv(path + f"output_tables/comm_{id}_{size}.tsv", sep="\t")
    for id in list_ids_analyzed:
        df = pd.read_csv(path + f"output_tables/comm_{id}_{size}.tsv", sep="\t")
        disease_name = dico_code_disease[str(id)]
        df = df.rename(columns={'Unnamed: 0': str(disease_name)})
        df = df.drop(['evidences'], axis=1)
        print(df)
        df.to_csv(path + f"output_tables/comm_{id}_{size}.tsv", sep="\t", index=False)

    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "output_tables/enrichment_communities_100.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

enrichment_communities(communities_100, list_ids_analyzed, 100)