import pandas as pd
import os
import sys

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../..')
os.chdir(path)
print(path)

from utilities import create_cluster_dico, filter_cluster

data_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data')
print(data_folder)
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
gmt_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data/gmt_files')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')

# set paths for GMT files
gmt_GOBP = os.path.join(gmt_folder, 'hsapiens.GO:BP.name.gmt')
print(gmt_GOBP)
gmt_GOCC = os.path.join(gmt_folder, 'hsapiens.GO:CC.name.gmt')
gmt_REAC = os.path.join(gmt_folder, 'hsapiens.REAC.name.gmt')

# create dictionary for clusters
dico_cluster_diseases = create_cluster_dico(cluster_output)
print(dico_cluster_diseases)
# filter the dico : only keep cluster containing at least 3 diseases
filtered_dico_cluster = filter_cluster(dico_cluster_diseases)
print(" ")
print(f"Clusters containing at least 3 diseases: {filtered_dico_cluster}")

def create_enrichment_files(cluster_id: str):
    """Function which creates enrichment files containing lists of 
    terms IDs for orsum analysis

    Args:
        cluster_id (str): the identifier of the cluster
    """
    os.mkdir(path + f'output_orsum/Orsum_{cluster_id}')
    dfEnrichmentGroup = pd.read_csv(path + f"output_tables/enrich_bioannot_{cluster_id}.tsv", sep="\t")
    sources=['GO:BP', 'GO:CC', 'REAC']
    for source in sources:
        dfEnrichmentGroupSource=dfEnrichmentGroup[dfEnrichmentGroup['source']==source]
        dfEnrichmentGroupSource=dfEnrichmentGroupSource.sort_values(by='p_value', ascending=True)
        dfEnrichmentGroupSource['native'].to_csv(path + f'output_orsum/Orsum_{cluster_id}/EnrichmentClust' + '-'+ source.replace(':','')+'.txt', index=False, header=None)


def applyOrsum(dico: dict, gmt: str, source: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
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
    # /!\ ORSUM PATH TO ADDAPT /!\
    # if orsum is installed as a conda pacgke
    #command = '/home/.../miniconda3/pkgs/orsum-1.6.0-hdfd78af_0/bin/orsum.py'
    command = '/home/cbeust/miniconda3/pkgs/orsum-1.6.0-hdfd78af_0/bin/orsum.py'
    # if orsum is installed in the current directory
    #command = path + 'orsum/orsum.py'
    command = command + ' --gmt \"'+gmt+'\" '
    command = command + '--files '
    for clusters in dico:
        filePath = summaryFolder+f'output_orsum/Orsum_{clusters}/EnrichmentClust-{source}.txt'
        print(f"FILEPATH : {filePath}")
        if os.stat(filePath).st_size == 0: # orsum exits when one enrichment file is empty
            noEnrichmentGroup.add(clusters)
        else:
            command=command+'\"'+filePath+'\" '
    command=command+'--fileAliases '
    for clusters in dico:
        if clusters not in noEnrichmentGroup:
            command=command+'\"'+str(clusters)+'\" '
    command=command+'--outputFolder \"'+outputFolder+'\" '
    command=command+'--maxRepSize '+ str(maxRepSize) + ' '
    command=command+'--minTermSize '+ str(minTermSize) + ' '
    command=command+'--numberOfTermsToPlot '+str(numberOfTermsToPlot)
    print('orsum command:')
    print(command)
    os.system(command)

def multiple_enrichment(dico_clusters: dict, source: str, gmt: str) -> None:
    os.mkdir(path + f"output_orsum/ME_results_{source}_clust")
    outputFolder = path + f"output_orsum/ME_results_{source}_clust"
    if source == "REAC":
        # for Reactome pathways we apply a different value for the maxRepSize parameter
        applyOrsum(dico=dico_clusters, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, maxRepSize=2000, numberOfTermsToPlot=50)
    else:
        applyOrsum(dico=dico_clusters, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, numberOfTermsToPlot=50)


def apply_orsum_to_clusters():
    for cluster in filtered_dico_cluster:
        create_enrichment_files(cluster)
    multiple_enrichment(dico_clusters=filtered_dico_cluster, source="GOBP", gmt=gmt_GOBP)
    multiple_enrichment(dico_clusters=filtered_dico_cluster, source="GOCC", gmt=gmt_GOCC)
    multiple_enrichment(dico_clusters=filtered_dico_cluster, source="REAC", gmt=gmt_REAC)

#apply_orsum_to_clusters()