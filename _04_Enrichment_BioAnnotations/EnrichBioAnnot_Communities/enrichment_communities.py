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
reac_genes = "../Build_background/REAC_genes.csv"
gobp_genes = "../Build_background/GOBP_genes.csv"
gocc_genes = "../Build_background/GOCC_genes.csv"

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

def compute_background_annot_mx(genes_file_annot: str, multiplex_nodes: list) -> list():
    df = pd.read_csv(genes_file_annot, sep=",", header=0)
    genes_annot = df["Gene"].to_list()
    background = set(genes_annot).intersection(set(multiplex_nodes))
    return list(background)

background_gobp = compute_background_annot_mx(genes_file_annot=gobp_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg go bp hugo: {len(background_gobp)}")
background_gocc = compute_background_annot_mx(genes_file_annot=gocc_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg go cc hugo: {len(background_gocc)}")
background_reac = compute_background_annot_mx(genes_file_annot=reac_genes, multiplex_nodes=all_nodes_HUGO)
print(f"bg reac hugo: {len(background_reac)}")

pa_diseases = pd.read_csv(orpha_names, sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

def enrichment_communities(source: list, list_comm: list, list_ids_analyzed: list, background: tuple) -> None:
    for comm, id in zip(list_comm, list_ids_analyzed):
        print(comm, id)
        with open(comm, 'r') as file:
            genes = []
            for line in file:
                genes += line.rsplit()
        gp = GProfiler(
        user_agent='https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/gost/profile/', 
        return_dataframe=True
        )
        enrich_gobp = gp.profile(
            organism='hsapiens', 
            query=genes,
            sources=[source[0]],
            domain_scope='custom',
            user_threshold=0.05,
            significance_threshold_method='false_discovery_rate', 
            background=background[0],
            no_evidences=False,
            )
        enrich_gocc = gp.profile(
            organism='hsapiens', 
            query=genes,
            sources=[source[1]],
            domain_scope='custom',
            user_threshold=0.05,
            significance_threshold_method='false_discovery_rate', 
            background=background[1],
            no_evidences=False,
            )
        enrich_reac = gp.profile(
            organism='hsapiens', 
            query=genes,
            sources=[source[2]],
            domain_scope='custom',
            user_threshold=0.05,
            significance_threshold_method='false_discovery_rate', 
            background=background[2],
            no_evidences=False,
            )
        enrich_gobp.to_csv(path + f"output_tables/comm_{id}_gobp.tsv", sep="\t")
        enrich_gocc.to_csv(path + f"output_tables/comm_{id}_gocc.tsv", sep="\t")
        enrich_reac.to_csv(path + f"output_tables/comm_{id}_reac.tsv", sep="\t")

enrichment_communities(source=["GO:BP", "GO:CC", "REAC"], list_comm=communities_100, list_ids_analyzed=list_ids_analyzed, background=(background_gobp, background_gocc, background_reac))

def merge_enrichment_files(list_ids_analyzed: list):
    for id in list_ids_analyzed:
        disease_name = dico_code_disease[str(id)]
        # read enrichment files
        enrich_gobp = pd.read_csv(path + f"output_tables/comm_{id}_gobp.tsv", sep="\t")
        enrich_gocc = pd.read_csv(path + f"output_tables/comm_{id}_gocc.tsv", sep="\t")
        enrich_reac = pd.read_csv(path + f"output_tables/comm_{id}_reac.tsv", sep="\t")
        # concatenate enrichment results in a single file
        all_annot = pd.concat([enrich_gobp, enrich_gocc, enrich_reac], axis=0)
        assert len(all_annot.index) == len(enrich_gobp.index) + len(enrich_gocc.index) + len(enrich_reac.index)
        # clean, sort and rename columns
        all_annot = all_annot.drop(['evidences'], axis=1)
        all_annot_sorted = all_annot.sort_values(by=["p_value"], ascending=True)
        all_annot_sorted.rename(columns={'p_value':'Corrected p_value', 'Unnamed: 0':f'{str(disease_name)}'}, inplace=True)
        all_annot_sorted.to_csv(f"output_tables/enrich_bioannot_comm_{id}.tsv", sep="\t", header=True, index=False)
        os.remove(f"output_tables/comm_{id}_gobp.tsv")
        os.remove(f"output_tables/comm_{id}_gocc.tsv")
        os.remove(f"output_tables/comm_{id}_reac.tsv")
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "output_tables/enrichment_bioannot_communities.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

merge_enrichment_files(list_ids_analyzed=list_ids_analyzed)