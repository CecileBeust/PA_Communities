# Import modules
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import argparse

def get_list_orpha_names(pa_diseases: str, id_diseases_analysed: list) -> list:
    """Function to get a list of PA diseases names
    from a file of diseases + their IDs and a list 
    of identifiers to analyze

    Args:
        path (str) : path of the working directory
        pa_diseases (str): _description_
        id_diseases_analysed (list): list of disease
        identifiers that we want to analyze

    Returns:
        list: list of disease names to analyze
    """
    disease_name_list = []
    file = pa_diseases
    with open(file, 'r') as fi:
        for line in fi:
            disease_id = line.split("\t")[0].rstrip()
            disease_name = line.split("\t")[1].rstrip()
            for id in id_diseases_analysed:        
                if id == disease_id[6:]:
                    disease_name_list.append(disease_name)
    return disease_name_list

def create_dico_disease_seeds(orpha_codes: str) -> tuple[dict, list]:
    """Function to create a dictionary of PA
    diseases and their associated causative
    genes (seeds)

    Args:
        path (str) : path of the working directory
        orpha_seeds (str): name of the file 
        containing the diseases IDs and their
        associated causative genes

    Returns:
        dict: dictionnary containing disease
        identifiers as keys and their seeds as
        values
        list : the list of disease identifiers
    """
    dico_seeds = {}
    list_disease = []
    file = orpha_codes
    with open(file, 'r') as fi:
        for line in fi:
            # separate diseases from seeds in the input file
            disease = line.split("\t")[0]
            list_disease.append(disease)
            # split seeds
            seeds_rsplit = line.split("\t")[1].rsplit()
            seeds = [genes.split(",") for genes in seeds_rsplit]
            # initialize key in dico for disease
            dico_seeds[disease] = []
            # writing one seeds file for each set of seeds
            # we take the orpha code of the disease to name the seeds files
            for list_genes in seeds:
                for genes in list_genes:
                    # add set of seeds in dico
                    dico_seeds[disease].append(genes)
        return (dico_seeds, list_disease)