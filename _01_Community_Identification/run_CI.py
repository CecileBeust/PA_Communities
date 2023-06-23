from itRWR import community_identification
import os

path = os.path.dirname(os.path.realpath(__file__))
path = path + '/'
os.chdir(path)

diseases = path + "orpha_codes_PA.txt"
num_iteration = 100
community_identification(path, diseases, num_iteration)