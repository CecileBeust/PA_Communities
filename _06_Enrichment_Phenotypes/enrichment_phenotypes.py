"""
Python script to perform an enrichment of clusters of communities
using HPO phenotypes (https://hpo.jax.org/app/)

We use an hypergeometric test to perform this enrichment. 
We consider the set of all ORPHANET diseases described in 
HPO as a background. For each HPO phenotype found associated 
to one or several disease(s) in a cluster, we compare the number 
of ORPHANET diseases associated to that phenotype in HPO, to the 
number of diseases in the cluster. 
"""

from gprofiler import GProfiler
import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import sys
from pyhpo import Ontology
_ = Ontology()

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../')
os.chdir(path)
print(path)

data_folder = os.path.join(os.path.dirname(__file__), '../', '_00_data')
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
clusters_compo = os.path.join(data_folder, 'clusters_100.tsv')

# variables statement
orpha_codes_file = pd.read_csv(orpha_codes, sep="\t", header=None)
orpha_names_file = pd.read_csv(orpha_names, sep="\t", header=None)

list_ids_analyzed = list()
pa_diseases_analyzed = list()
all_pa_diseases = list()
dico_code_disease = dict()

# extract diseases ORPHANET codes analyzed
for index, row in orpha_codes_file.iterrows():
    if not row[1] is np.nan and not pd.isnull(row[1]):
        list_ids_analyzed.append(row[0])

print(list_ids_analyzed)
print(f"{len(list_ids_analyzed)} diseases analyzed ")
print(" ")

# create list of disease names
for index, row in orpha_names_file.iterrows():
    all_pa_diseases.append(row[1])
    if int(row[0][6:]) in list_ids_analyzed:
        pa_diseases_analyzed.append(row[1])

print(f"{len(pa_diseases_analyzed)} PA diseases: {pa_diseases_analyzed}")
print(" ")

# create dico of ORPHANET codes and diseases names
for index, row in orpha_names_file.iterrows():
    code = row[0][6:]
    disease = row[1]
    dico_code_disease[code] = disease

print("Dictionary ORPHANET diseases")
print(dico_code_disease)
print(" ")

def create_dico_clusters_diseases(clusters_file) -> dict:
    """Function that creates a dictionary of clusters
    and the diseases they contain

    Args:
        clusters_file (str): name of the file where
        the clusters are detailed

    Returns:
        dict: dico where keys are clusters and values
        are the diseases they contain
    """
    df = pd.read_csv(clusters_file, sep="\t", header=0)
    dico_clusters_diseases = dict()
    for index, row in df.iterrows():
        if not row[0] in dico_clusters_diseases.keys():
            dico_clusters_diseases[row[0]] = [row[1]]
        else:
            dico_clusters_diseases[row[0]] += [row[1]]
    return dico_clusters_diseases
            
dico_clusters_diseases = create_dico_clusters_diseases(clusters_compo)
print(f"Dico clusters diseases : {dico_clusters_diseases}")
print(" ")

cluster1 = dico_clusters_diseases[1]
cluster2 = dico_clusters_diseases[2]
cluster3 = dico_clusters_diseases[3]
cluster4 = dico_clusters_diseases[4]
cluster5 = dico_clusters_diseases[5]
cluster6 = dico_clusters_diseases[6]
list_clusters = [cluster1, cluster2, cluster3, cluster4, cluster5, cluster6]
print(list_clusters)

def build_background(list_diseases: list, dico_code_disease: dict) -> set:
    """Function that determines the background of HPO phenotypes to use
    for the hypergeometric test

    Args:
        list_diseases (list): PA diseases
        dico_diseases_code (dict): dico of PA diseases
        and their associated ORPHANET codes

    Returns:
        list: set containing all the phenotypes of the background
    """
    background = list()
    for disease in list_diseases:
        disease_id = {id for id in dico_code_disease if dico_code_disease[id]==disease}
        for id in disease_id:
            pheno_file = path + f"Data_HPO/terms_for_ORPHA_{id}.xlsx"
            df = pd.read_excel(pheno_file, header=0)
            phenotypes = df['HPO_TERM_NAME'].to_list()
            background.extend(phenotypes)
    return set(background)

background = build_background(pa_diseases_analyzed, dico_code_disease)
print(" ")
print(f"Background : {len(background)} phenotypes")

def create_dico_clusters_pheno(list_diseases_in_clusters: list, dico_code_disease: dict) -> dict:
    """Function which creates a dictionary of HPO phenotypes associated
    to PA diseases

    Args:
        list_diseases_in_clusters (list): list of diseases in clusters

    Returns:
        dict: dico of clusters as keys and their associated phenotypes 
        from HPO as values
    """
    dico_clusters_pheno = dict()
    i = 1
    for cluster in list_diseases_in_clusters:
        # initialize list of phenotypes of the cluster
        all_phenotypes = []
        for disease in cluster:
            disease_id = {id for id in dico_code_disease if dico_code_disease[id]==disease}
            print(disease_id)
            for id in disease_id:
                pheno_file =  path + f"Data_HPO/terms_for_ORPHA_{id}.xlsx"
                print(pheno_file)
            df = pd.read_excel(pheno_file, header=0)
            phenotypes = df['HPO_TERM_NAME'].to_list()
            all_phenotypes.extend(phenotypes)
            print(all_phenotypes)
        dico_clusters_pheno[i] = all_phenotypes
        i += 1
    return dico_clusters_pheno

"""dico_clusters_pheno = create_dico_clusters_pheno(list_clusters, dico_code_disease)
print(" ")
print(dico_clusters_pheno)"""

def enrich_phenotypes(background_pheno: set, diseases_pa: list, list_cluster: list, background_diseases_HPO: int):
    j = 1
    for cluster in list_cluster:
        df = pd.DataFrame(columns=['Cluster', 'HPO ID', 'Phenotype', 'p-value', 'Corrected p-value'])
        nb_diseases_cluster = len(cluster)
        i = 0
        for phenotype in background_pheno:
            term = Ontology.get_hpo_object(phenotype)
            list_diseases_pheno = list()
            for disease in term.orpha_diseases:
                #if str(disease) in diseases_pa:
                list_diseases_pheno.append(str(disease))
            nb_diseases_pheno = len(set(list_diseases_pheno))
            success_sample = len(set(cluster).intersection(set(list_diseases_pheno)))
            p_value = stats.hypergeom.sf(success_sample-1, background_diseases_HPO, nb_diseases_pheno, nb_diseases_cluster)
            df._set_value(i, 'Cluster', j)
            df._set_value(i, 'HPO ID', str(term)[:10])
            df._set_value(i, 'Phenotype', phenotype)
            df._set_value(i, 'p-value', p_value)
            i += 1
        pvals = df['p-value'].to_list()
        p_adjusted = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        df['Corrected p-value'] = p_adjusted
        df.to_csv(path + f"output_tables/enrichment_cluster{j}_all_pheno.tsv", sep="\t", index=False)
        j += 1
    tsv_dir = Path(path + "output_tables/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "output_tables/enrichment_all_phenotypes.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

enrich_phenotypes(background_pheno = background, diseases_pa = all_pa_diseases, list_cluster = list_clusters, background_diseases_HPO = 4262)

def create_gmt():
    onto = pd.read_csv(path + "Data_HPO/HP.csv", sep=",", header=0)
    print(onto)
    to_write = list()
    for index, row in onto.iterrows():
        line = str()
        if row[4] == True:
            if not row[54] is np.nan or not pd.isnull(row[54]):
                id_hpo = str(row[54][31:]).replace("_",":")
        else:
            id_hpo = str(row[0][31:]).replace("_",":")
            name = row[1]
            term = Ontology.get_hpo_object(id_hpo)
            diseases = term.orpha_diseases
            line += (str(id_hpo) + "\t" + str(name) + "\t")
            for dis in diseases:
                line += (str(dis) + "\t")
            to_write += [line]
    with open(path + "hpo_phenotypes.gmt", 'w') as file:
        file.write("\n".join([line for line in to_write]))

create_gmt()

def create_enrichment_files():
    for i in range(1,7):
        os.mkdir(path + f"output_orsum/Orsum_cluster_{i}")
        dfEnrichmentGroup = pd.read_csv(path + f"output_tables/enrichment_cluster{i}_all_pheno.tsv", sep="\t", header=0)
        dfEnrichmentGroup = dfEnrichmentGroup.sort_values(by="Corrected p-value", ascending=True)
        dfEnrichmentGroup = dfEnrichmentGroup[dfEnrichmentGroup["Corrected p-value"] <= 0.05]
        print(dfEnrichmentGroup)
        dfEnrichmentGroupSource = dfEnrichmentGroup["HPO ID"]
        print(dfEnrichmentGroupSource)
        dfEnrichmentGroupSource.to_csv(path + f"output_orsum/Orsum_cluster_{i}/Enrich-HPO-cluster-{i}.txt", index=False, header=None)

create_enrichment_files()

def applyOrsum(dico: dict, gmt: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
    """Function to apply orsum on enrichment analysis results of clusters inside a dictionnary

    Args:
        dico (dict): the dictionnary containing the clusters ID and the associated diseases
        size (int): size of communities analyzed (30, 50, 100)
        th (float): threshold used to determine the clusters
        gmt (str): name of the GMT file used
        source (str): enrichment result source analyzed
        summaryFolder (str): folder name for storing the results of the analysis
        outputFolder (str): folder name for storing orsum outputs
        maxRepSize (_type_, optional): maximal size of a representative term. Defaults to int(1E6).
        minTermSize (int, optional): The minimum size of the terms to be processed. Defaults to 10.
        numberOfTermsToPlot (int, optional): The number of representative terms to be presented in barplot and heatmap. Defaults to 50.
    """
    noEnrichmentGroup = set()
    #/!\ ORSUM PATH TO ADDAPT /!\
    # if orsum is installed as a conda pacgke
    command = '/home/cbeust/miniconda3/pkgs/orsum-1.6.0-hdfd78af_0/bin/orsum.py'
    # if orsum is installed in the current directory
    #command = path + 'orsum/orsum.py'
    command = command + ' --gmt \"'+gmt+'\" '
    command = command + '--files '
    for clusters in dico:
        filePath = summaryFolder+f'output_orsum/Orsum_cluster_{int(clusters)}/Enrich-HPO-cluster-{int(clusters)}.txt'
        if os.stat(filePath).st_size == 0: # orsum exits when one enrichment file is empty
            noEnrichmentGroup.add(int(clusters))
        else:
            command=command+'\"'+filePath+'\" '
    command=command+'--fileAliases '
    for clusters in dico:
        if int(clusters) not in noEnrichmentGroup:
            command=command+'\"'+str(int(clusters))+'\" '
    command=command+'--outputFolder \"'+outputFolder+'\" '
    command=command+'--maxRepSize '+ str(maxRepSize) + ' '
    command=command+'--minTermSize '+ str(minTermSize) + ' '
    command=command+'--numberOfTermsToPlot '+str(numberOfTermsToPlot)
    print('orsum command:')
    print(command)
    os.system(command)


def multiple_enrichment(dico_clusters_diseases: dict, gmt: str):
    os.mkdir(path + f"output_orsum/ME_results_phenotypes")
    outputFolder = path + f"output_orsum/ME_results_phenotypes"
    applyOrsum(dico=dico_clusters_diseases, gmt=gmt, summaryFolder=path, outputFolder=outputFolder)


# set paths for GMT files
gmt_HPO = path + "hpo_phenotypes.gmt"

multiple_enrichment(dico_clusters_diseases, gmt_HPO)

