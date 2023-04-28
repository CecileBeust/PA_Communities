import pandas as pd
import os
from enrichment_communities import *
import collections

path = "/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/Analysis_Communities_V2/enrichment"


def create_cluster_dico(cluster_file: str):
    """Function to create a dictionary containing clusters as keys
    and the associated disease list as values

    Args:
        cluster_file (str): the name of the file containing cluster assignements

    Returns:
        dict: dictionnary of cluster-diseases associations
    """
    df = pd.read_csv(cluster_file, sep="\t")
    dico_cluster_diseases = {}
    i = 0
    for cluster in df['cluster']:
        disease = df.iloc[i]['disease']
        if cluster not in dico_cluster_diseases.keys():
            dico_cluster_diseases[cluster] = [disease]
        else:
            dico_cluster_diseases[cluster] += [disease]
        i += 1
    return dico_cluster_diseases


def filter_cluster(dico_cluster: dict):
    """Function to filter a dictionary of cluster assignments
    to only keep clusters containing at least 3 diseases

    Args:
        dico_cluster (dict): the dictionary obtained with the 
        create_cluster_dico() function

    Returns:
        dict: a filtered dictionary
    """
    filtered_dict = {}
    for cluster in dico_cluster:
        if len(dico_cluster[cluster]) >=3 :
            filtered_dict[cluster] = dico_cluster[cluster]
    return filtered_dict


def create_enrichment_files(size: int, cluster_id: int, th: float):
    """Function which creates enrichment files containing lists of 
    terms IDs for orsum analysis

    Args:
        size (int): the size of communities analyzed (30, 50, 100)
        cluster_id (int): the ID of the cluster
        th (float): threshold used to determine the clusters
    """
    os.mkdir(path + f"/{size}_{th}/" + f"Orsum_{size}_{th}_cluster_{cluster_id}")
    dfEnrichmentGroup = pd.read_csv(path + "/" + f"{size}_{th}/cluster_{size}_{cluster_id}.tsv", sep="\t")
    sources=['GO:BP', 'GO:CC', 'REAC']
    for source in sources:
        dfEnrichmentGroupSource=dfEnrichmentGroup[dfEnrichmentGroup['source']==source]
        dfEnrichmentGroupSource=dfEnrichmentGroupSource.sort_values(by='p_value', ascending=True)
        dfEnrichmentGroupSource['native'].to_csv(path + f'/{size}_{th}/' + f'Orsum_{size}_{th}_cluster_{cluster_id}/EnrichmentComm' + '-'+ source.replace(':','')+'.txt', index=False, header=None)


def applyOrsum2(dico: dict, size: int, th: float, gmt: str, source: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
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
    #command = '/home/cbeust/miniconda3/pkgs/orsum-1.4.0-hdfd78af_0/bin/orsum.py'
    command = '/home/cbeust/Landscape_PA/CommunityIdentification/CommunityIdentification_V2/Analysis_Communities_V2/orsum/orsum.py'
    command = command + ' --gmt \"'+gmt+'\" '
    command = command + '--files '
    for clusters in dico:
        filePath = summaryFolder+os.sep+f'{size}_{th}/Orsum_{size}_{th}_cluster_{int(clusters)}/EnrichmentComm-{source}.txt'
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


# set paths for GMT files
gmt_GOBP = path + "/" + "gmt_files/hsapiens.GO:BP.name.gmt"
gmt_GOCC = path + "/" + "gmt_files/hsapiens.GO:CC.name.gmt"
gmt_REAC = path + "/" + "gmt_files/hsapiens.REAC.name.gmt" 

# create dictionary for clusters
dico_cluster_diseases_100_0_7 = create_cluster_dico("cluster_output_100_0.7.tsv")
print(dico_cluster_diseases_100_0_7)
# filter the dico : only keep cluster containing at least 3 diseases
filtered_dico_cluster_100_0_7 = filter_cluster(dico_cluster_diseases_100_0_7)
print(filtered_dico_cluster_100_0_7)
# sort the dico to have the right order on the orsum heatmap
sorted_dico_clusters_100_0_7 = collections.OrderedDict(sorted(filtered_dico_cluster_100_0_7.items()))

# create enrichment files
"""for cluster in sorted_dico_clusters_100_0_7:
    create_enrichment_files(100, int(cluster), 0.7)"""

def multiple_enrichment(th, size, source, gmt):
    os.mkdir(path + f"/{size}_{th}/ME_results_{source}_comm_{th}")
    outputFolder = path + f"/{size}_{th}/ME_results_{source}_comm_{th}"
    applyOrsum2(dico=sorted_dico_clusters_100_0_7, size=size, th=0.7, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, maxRepSize=2000, numberOfTermsToPlot=50)


"""multiple_enrichment(0.7, size=100, source="GOBP", gmt=gmt_GOBP)
multiple_enrichment(0.7, size=100, source="GOCC", gmt=gmt_GOCC)
multiple_enrichment(0.7, size=100, source="REAC", gmt=gmt_REAC)"""