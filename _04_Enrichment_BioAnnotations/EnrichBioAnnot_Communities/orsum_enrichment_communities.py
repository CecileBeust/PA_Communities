import pandas as pd
import os
import sys
import argparse

# define path
path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
sys.path.append('../..')
os.chdir(path)
print(path)

from utilities import create_dico_disease_seeds, build_communities_list

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

data_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data')
print(data_folder)
orpha_codes = os.path.join(data_folder, 'orpha_codes_PA.txt')
orpha_names = os.path.join(data_folder, 'pa_orphanet_diseases.tsv')
gmt_folder = os.path.join(os.path.dirname(__file__), '../..', '_00_data/gmt_files')
cluster_output = os.path.join(data_folder, 'cluster_output_100_0.7.tsv')

# set paths for GMT files
gmt_GOBP = os.path.join(gmt_folder, 'hsapiens.GO:BP.name.gmt')
gmt_GOCC = os.path.join(gmt_folder, 'hsapiens.GO:CC.name.gmt')
gmt_REAC = os.path.join(gmt_folder, 'hsapiens.REAC.name.gmt')


(dico_disease_seeds, list_id) = create_dico_disease_seeds(orpha_codes)
(communities_100, not_analyzed) = build_communities_list(comm_path, list_id, 100)
list_ids_analyzed = [x for x in list_id if x not in not_analyzed]

def create_enrichment_files(comm_id: str):
    """Function which creates enrichment files containing lists of 
    terms IDs for orsum analysis

    Args:
        cluster_id (str): the identifier of the cluster
    """
    os.mkdir(path + f'output_orsum/Orsum_{comm_id}')
    dfEnrichmentGroup = pd.read_csv(path + f"output_tables/comm_{comm_id}_100.tsv", sep="\t")
    sources=['GO:BP', 'GO:CC', 'REAC']
    for source in sources:
        dfEnrichmentGroupSource=dfEnrichmentGroup[dfEnrichmentGroup['source']==source]
        dfEnrichmentGroupSource=dfEnrichmentGroupSource.sort_values(by='p_value', ascending=True)
        dfEnrichmentGroupSource['native'].to_csv(path + f'output_orsum/Orsum_{comm_id}/EnrichmentClust' + '-'+ source.replace(':','')+'.txt', index=False, header=None)


def applyOrsum(list_comm: list, gmt: str, source: str, summaryFolder: str, outputFolder: str, maxRepSize=int(1E6), minTermSize=10, numberOfTermsToPlot=50):
    """Function to apply orsum on enrichment analysis results of clusters inside a dictionnary

    Args:
        dico (dict): the dictionnary containing the clusters ID and the associated diseases
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
    for comm in list_comm:
        print(comm)
        filePath = summaryFolder+f'output_orsum/Orsum_{comm}/EnrichmentClust-{source}.txt'
        print(filePath)
        if os.stat(filePath).st_size == 0: # orsum exits when one enrichment file is empty
            noEnrichmentGroup.add(comm)
        else:
            command=command+'\"'+filePath+'\" '
    command=command+'--fileAliases '
    for comm in list_comm:
        if comm not in noEnrichmentGroup:
            command=command+'\"'+str(comm)+'\" '
    command=command+'--outputFolder \"'+outputFolder+'\" '
    command=command+'--maxRepSize '+ str(maxRepSize) + ' '
    command=command+'--minTermSize '+ str(minTermSize) + ' '
    command=command+'--numberOfTermsToPlot '+str(numberOfTermsToPlot)
    print('orsum command:')
    print(command)
    os.system(command)

def multiple_enrichment(list_comm: list, source: str, gmt: str) -> None:
    os.mkdir(path + f"output_orsum/ME_results_{source}_comm")
    outputFolder = path + f"output_orsum/ME_results_{source}_comm"
    if source == "REAC":
        # for Reactome pathways we apply a different value for the maxRepSize parameter
        applyOrsum(list_comm=list_ids_analyzed, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, maxRepSize=2000, numberOfTermsToPlot=50)
    else:
        applyOrsum(list_comm=list_ids_analyzed, gmt=gmt, source=source, summaryFolder=path, outputFolder=outputFolder, numberOfTermsToPlot=50)


def apply_orsum_to_clusters():
    for comm in list_ids_analyzed:
        create_enrichment_files(comm)
    multiple_enrichment(list_comm=list_ids_analyzed, source="GOBP", gmt=gmt_GOBP)
    multiple_enrichment(list_comm=list_ids_analyzed, source="GOCC", gmt=gmt_GOCC)
    multiple_enrichment(list_comm=list_ids_analyzed, source="REAC", gmt=gmt_REAC)

apply_orsum_to_clusters()