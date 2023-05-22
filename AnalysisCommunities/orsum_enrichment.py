import pandas as pd
import os
from enrichment_communities import create_cluster_dico, filter_cluster
import collections

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)
print(path)


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


def applyOrsum(dico: dict, size: int, th: float, gmt: str, source: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
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
    # if orsum is installed as a conda pacgke
    #command = '/home/.../miniconda3/pkgs/orsum-1.4.0-hdfd78af_0/bin/orsum.py'
    # if orsum is installed in the current directory
    #command = path + 'orsum/orsum.py'
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
gmt_GOBP = path + "gmt_files/hsapiens.GO:BP.name.gmt"
gmt_GOCC = path + "gmt_files/hsapiens.GO:CC.name.gmt"
gmt_REAC = path + "gmt_files/hsapiens.REAC.name.gmt" 

# create dictionary for clusters
dico_cluster_diseases = create_cluster_dico(path + "cluster_output_100_0.7.tsv")
print(dico_cluster_diseases)
# filter the dico : only keep cluster containing at least 3 diseases
filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(filtered_dico_cluster)
# sort the dico to have the right order on the orsum heatmap
sorted_dico_clusters = collections.OrderedDict(sorted(filtered_dico_cluster.items()))

def create_enrichment_files_clusters(sorted_dico_clusters: dict, size: int) -> None:
    for cluster in sorted_dico_clusters:
        create_enrichment_files(size, int(cluster), 0.7)

def multiple_enrichment(sorted_dico_clusters: dict, th: float, size: int, source: str, gmt: str) -> None:
    os.mkdir(path + f"/{size}_{th}/ME_results_{source}_comm_{th}")
    outputFolder = path + f"/{size}_{th}/ME_results_{source}_comm_{th}"
    if source == "REAC":
        applyOrsum(dico=sorted_dico_clusters, size=size, th=0.7, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, maxRepSize=2000, numberOfTermsToPlot=50)
    else:
        applyOrsum(dico=sorted_dico_clusters, size=size, th=0.7, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, numberOfTermsToPlot=50)


def apply_orsum_to_clusters(th: float, size: int):
    multiple_enrichment(th=th, size=size, source="GOBP", gmt=gmt_GOBP)
    multiple_enrichment(th=th, size=size, source="GOCC", gmt=gmt_GOCC)
    multiple_enrichment(th=th, size=size, source="REAC", gmt=gmt_REAC)

apply_orsum_to_clusters(0.7)