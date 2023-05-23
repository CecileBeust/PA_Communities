import pandas as pd
from analyse_final_communities import create_dico_disease_seeds
from cluster_communities import build_communities_list
from enrichment_clusters import create_cluster_dico
from gprofiler import GProfiler
import os
import ast
import argparse
from pathlib import Path

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

# Argparse
parser = argparse.ArgumentParser(
    prog="cluster_communities.py", 
    description="functions cluster communities based on Jaccard index"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

# Check if the path exist
if os.path.exists(comm_path) == False :
    raise ValueError("Incorrect path, please try again")


(dico_disease_seeds, list_id) = create_dico_disease_seeds(path, "data/orpha_codes_PA.txt")
(communities_10, not_analyzed) = build_communities_list(comm_path, list_id, 10)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]

pa_diseases = pd.read_csv(path + 'data/pa_orphanet_diseases.tsv', sep="\t", header=None)
dico_code_disease = {}
for index, row in pa_diseases.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease
print(dico_code_disease)

def enrichment_communities(list_comm, list_ids_analyzed, size):
    print(len(list_comm))
    for comm, id in zip(list_comm, list_ids_analyzed):
        with open(comm, 'r') as file:
            genes = []
            for line in file:
                genes += line.rsplit()
            print(genes)
            gp = GProfiler(return_dataframe=True)
            enrich = gp.profile(organism='hsapiens', query=genes, no_evidences=False)
            print(enrich)
            enrich.to_csv(f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3/Analysis_Communities_V3/enrichment/EnrichCommunities/comm_{id}_{size}.tsv", sep="\t")
    for id in list_ids_analyzed:
        df = pd.read_csv(f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3/Analysis_Communities_V3/enrichment/EnrichCommunities/comm_{id}_{100}.tsv", sep="\t")
        disease_name = dico_code_disease[str(id)]
        df = df.rename(columns={'Unnamed: 0': str(disease_name)})
        df = df.drop(['evidences'], axis=1)
        print(df)
        df.to_csv(f"/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3/Analysis_Communities_V3/enrichment/EnrichCommunities/comm_{id}_{100}.tsv", sep="\t", index=False)

    tsv_dir = Path("/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3/Analysis_Communities_V3/enrichment/EnrichCommunities/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter('/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V3/Analysis_Communities_V3/enrichment/EnrichCommunities/enrichment_communities_100.xlsx', engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

enrichment_communities(communities_10, list_ids_analyzed, 10)