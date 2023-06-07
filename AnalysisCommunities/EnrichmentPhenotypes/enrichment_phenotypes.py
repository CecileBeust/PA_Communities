from gprofiler import GProfiler
import pandas as pd
import numpy as np
import scipy.stats as stats
from pathlib import Path
from statsmodels.sandbox.stats.multicomp import multipletests
import os
import seaborn as sns
import matplotlib.pyplot as plt
from pyhpo import Ontology
import argparse
_ = Ontology()

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)

diseases_pa = []
dico_diseases_code = dict()
pa = pd.read_csv("PA_diseases.tsv", sep="\t", header=0)
print(pa)
for index, row in pa.iterrows():
    if not row[2] is np.nan or not pd.isnull(row[2]):
        diseases_pa.append(row[0])
        dico_diseases_code[row[0]] = row[1]

# remove this disease because seed not in network
diseases_pa.remove("Arterial tortuosity syndrome")
print(f"PA diseases : {diseases_pa}")
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
            
dico_clusters_diseases = create_dico_clusters_diseases("clusters_100.tsv")
print(f"Dico clusters diseases : {dico_clusters_diseases}")
print(" ")

def hypergeome(list1, list2, background) -> float:
    """Function to perform hypergeometric test between
    two lists of terms

    Args:
        list1 (list): first list of terms
        list2 (list): second list of terms
        background (int): number of terms taken as a 
        background for the statistical significance
        of the analysis

    Returns:
        floast: p-value of the hypergeometric test
    """
    # Define the number of features in the two lists
    list1_size = len(list1)
    list2_size = len(list2)

    # Determine the number of genes that are common to both lists
    common = set(list1).intersection(list2)
    common_size = len(common)

    # Define the number of features to randomly sample from the background
    sample_size = list1_size

    # Calculate the p-value using a hypergeometric test
    p_value = stats.hypergeom.sf(common_size-1, background, list2_size, sample_size)

    #print(f'p-value: {p_value:.4f}')
    return p_value


def get_PA_pheno_for_disease(clusters_file: str) -> tuple:
    """Extract diseases from each cluster

    Args:
        clusters_file (str): cluster assignement file

    Returns:
        tuple: tuple of list of diseases in each cluster
    """
    clusters = pd.read_csv(clusters_file, sep="\t", header=0)
    cluster1 = []
    cluster2 = []
    cluster3 = []
    cluster4 = []
    cluster5 = []
    cluster6 = []
    for i in range(1,7):
        for index, row in clusters.iterrows():
            if row[0] == i:
                if i == 1:
                    cluster1.append(row[1])
                elif i == 2:
                    cluster2.append(row[1])
                elif i == 3:
                    cluster3.append(row[1])
                elif i == 4:
                    cluster4.append(row[1])
                elif i == 5:
                    cluster5.append(row[1])
                elif i == 6:
                    cluster6.append(row[1])
    return (cluster1, cluster2, cluster3, cluster4, cluster5, cluster6)
                
(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6) = get_PA_pheno_for_disease("clusters_100.tsv")
list_clusters = [cluster1, cluster2, cluster3, cluster4, cluster5, cluster6]
print(len(list_clusters))

def enrich__PA_phenotypes(list_clusters: list, dico_pheno: dict, background: int):
    df_pval = pd.DataFrame(columns=['Cluster', 'DisCluster', 'Phenotype', 'Diseases', 'Intersection', 'pvalue', 'Corrected pvalue'])
    j = 1
    for cluster in list_clusters:
        print(j)
        i = 0
        for pheno, diseases in dico_pheno.items():
            p = hypergeome(cluster, diseases, background)
            df_pval._set_value(i, 'Cluster', j)
            df_pval._set_value(i, 'DisCluster', cluster)
            df_pval._set_value(i, 'Phenotype', pheno)
            df_pval._set_value(i, 'Diseases', diseases)
            df_pval._set_value(i, 'Intersection', set(cluster).intersection(diseases))
            df_pval._set_value(i, 'pvalue', p)
            i += 1
            pvals = df_pval['pvalue'].to_list()
        p_adjusted = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        df_pval['Corrected pvalue'] = p_adjusted
        df_pval.to_csv(path + f"/enrichment_cluster{j}_phenotypes.tsv", sep="\t", index=False)  
        j += 1
    
    tsv_dir = Path(path)
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "/enrichment_phenotypes.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

#enrich__PA_phenotypes(list_clusters, dico_pheno_diseases, 67)

def build_background(list_diseases: list, dico_diseases_code: dict) -> list:
    """Function that determines the background of HPO phenotypes to use
    for the hypergeometric test

    Args:
        list_diseases (list): PA diseases
        dico_diseases_code (dict): dico of PA diseases
        and their associated ORPHANET codes

    Returns:
        list: list containing all the phenotypes of the background
    """
    background= list()
    for disease in list_diseases:
        disease_id = dico_diseases_code[disease]
        pheno_file = path + f"HPOTermsOrphaDiseases/terms_for_ORPHA_{disease_id}.xlsx"
        df = pd.read_excel(pheno_file, header=0)
        phenotypes = df['HPO_TERM_NAME'].to_list()
        background.extend(phenotypes)
    return background

background = build_background(diseases_pa, dico_diseases_code)
print(f"Background : {len(background)} phenotypes")

def create_dico_clusters_pheno(list_diseases_in_clusters: list) -> dict:
    """Function which creates a dictionary of HPO phenotypes associated
    to PA diseases

    Args:
        list_diseases_in_clusters (list): list of diseases in clusters

    Returns:
        dict: dico of diseases as keys and their associated phenotypes 
        from HPO as values
    """
    dico_clusters_pheno = dict()
    i = 1
    for cluster in list_diseases_in_clusters:
        all_phenotypes = []
        for disease in cluster:
            disease_id = dico_diseases_code[disease]
            pheno_file =  path + '/' + f"HPOTermsOrphaDiseases/terms_for_ORPHA_{disease_id}.xlsx"
            df = pd.read_excel(pheno_file, header=0)
            phenotypes = df['HPO_TERM_NAME'].to_list()
            all_phenotypes.extend(phenotypes)
        dico_clusters_pheno[i] = all_phenotypes
        i += 1
    return dico_clusters_pheno

dico_clusters_pheno = create_dico_clusters_pheno(list_clusters)
print(" ")
print(len(dico_clusters_pheno))

def enrich_phenotypes(background_pheno: list, diseases_pa: list, list_cluster: list, back: int):
    j = 1
    for cluster in list_cluster:
        df = pd.DataFrame(columns=['Cluster', 'HPO ID', 'Phenotype', 'p-value', 'Corrected p-value'])
        nb_diseases_cluster = len(cluster)
        i = 0
        for phenotype in set(background_pheno):
            term = Ontology.get_hpo_object(phenotype)
            list_diseases_pheno = list()
            for disease in term.orpha_diseases:
                #if str(disease) in diseases_pa:
                list_diseases_pheno.append(str(disease))
            nb_diseases_pheno = len(set(list_diseases_pheno))
            success_sample = len(set(cluster).intersection(set(list_diseases_pheno)))
            p_value = stats.hypergeom.sf(success_sample-1, back, nb_diseases_pheno, nb_diseases_cluster)
            df._set_value(i, 'Cluster', j)
            df._set_value(i, 'HPO ID', str(term)[:10])
            df._set_value(i, 'Phenotype', phenotype)
            df._set_value(i, 'p-value', p_value)
            i += 1
        pvals = df['p-value'].to_list()
        p_adjusted = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        df['Corrected p-value'] = p_adjusted
        df.to_csv(path + f"EnrichPhenoResults/enrichment_cluster{j}_all_pheno.tsv", sep="\t", index=False)
        j += 1
    tsv_dir = Path(path + "EnrichPhenoResults/")
    tsv_data = {}
    for tsv_file in tsv_dir.glob('*.tsv'):
        tsv_name = tsv_file.stem
        tsv_data[tsv_name] = pd.read_csv(tsv_file, sep="\t")
    writer = pd.ExcelWriter(path + "EnrichPhenoResults/enrichment_all_phenotypes.xlsx", engine='xlsxwriter')
    for sheet_name, sheet_data in tsv_data.items():
        sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)
    writer.save()

pa = pd.read_csv("PA_diseases.tsv", sep="\t", header=None)
diseases_pa_complete = pa[0].to_list()

enrich_phenotypes(background_pheno = background, diseases_pa = diseases_pa_complete, list_cluster = list_clusters, back = 4262)

def create_gmt():
    onto = pd.read_csv(path + "HP.csv", sep=",", header=0)
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
        os.mkdir(path + f"Orsum_cluster_{i}")
        dfEnrichmentGroup = pd.read_csv(path + f"EnrichPhenoResults/enrichment_cluster{i}_all_pheno.tsv", sep="\t", header=0)
        dfEnrichmentGroup = dfEnrichmentGroup.sort_values(by="Corrected p-value", ascending=True)
        dfEnrichmentGroup = dfEnrichmentGroup[dfEnrichmentGroup["Corrected p-value"] <= 0.05]
        print(dfEnrichmentGroup)
        dfEnrichmentGroupSource = dfEnrichmentGroup["HPO ID"]
        print(dfEnrichmentGroupSource)
        dfEnrichmentGroupSource.to_csv(path + f"/Orsum_cluster_{i}/Enrich-HPO-cluster-{i}.txt", index=False, header=None)

create_enrichment_files()

def applyOrsum2(dico: dict, gmt: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
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
    command = '/home/cbeust/miniconda3/pkgs/orsum-1.6.0-hdfd78af_0/bin/orsum.py'
    #command = 'orsum.py'
    command = command + ' --gmt \"'+gmt+'\" '
    command = command + '--files '
    for clusters in dico:
        filePath = summaryFolder+os.sep+f'/Orsum_cluster_{int(clusters)}/Enrich-HPO-cluster-{int(clusters)}.txt'
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


def multiple_enrichment(dico_clusters_diseases, gmt):
    os.mkdir(path + f"/ME_results_phenotypes")
    outputFolder = path + f"/ME_results_phenotypes"
    applyOrsum2(dico=dico_clusters_diseases, gmt=gmt, summaryFolder=path, outputFolder=outputFolder)


# set paths for GMT files
gmt_HPO = path + "/hpo_phenotypes.gmt"

multiple_enrichment(dico_clusters_diseases, gmt_HPO)

