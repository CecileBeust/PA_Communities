import os
import sys
import argparse
import pandas as pd

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../..')
os.chdir(path)
print(path)

from utilities import load_networks
reac_genes = "REAC_genes.csv"
gobp_genes = "GOBP_genes.csv"
gocc_genes = "GOCC_genes.csv"

# Argparse
parser = argparse.ArgumentParser(
    prog="enrichment_clusters.py", 
    description="functions perform enrichment of clusters of gene communities"
    )
parser.add_argument("-p", "--path", help="path where communities are stored", required=True, type=str)
parser.add_argument("-v", "--verbose", help="increase output verbosity", required=False, action="store_true")
args = parser.parse_args()

comm_path = args.path

all_nodes_HUGO = list(load_networks(comm_path=comm_path))

def compute_background_annot_mx(genes_file_annot: str, multiplex_nodes: list) -> list():
    df = pd.read_csv(genes_file_annot, sep=",", header=0)
    print(df)
    genes_annot = df["Gene"].to_list()
    background = set(genes_annot).intersection(set(multiplex_nodes))
    print(len(background))
    return background

background_gobp = compute_background_annot_mx(genes_file_annot=gobp_genes, multiplex_nodes=all_nodes_HUGO)
background_go_cc = compute_background_annot_mx(genes_file_annot=gocc_genes, multiplex_nodes=all_nodes_HUGO)
background_reac = compute_background_annot_mx(genes_file_annot=reac_genes, multiplex_nodes=all_nodes_HUGO)