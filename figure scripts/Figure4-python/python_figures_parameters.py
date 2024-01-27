# -*- coding: utf-8 -*-
"""

Parameters used in the creation of the low-dimensional analysis subplots in Figure 4, 
and the Supporting Information figure 11 for the paper "Efficient encoding of aversive location by CA3 long-range projections"

Created by Albert Miguel López

                                      
"""


import os.path

##### Paths #####
#Add Path to the general folder for the project
recognized_project_paths = [
    "D:\\MPI_Brain\\Other Projects\\Negar Paper\\Data for python",
    "D:\Albert\MPI_Brain\Bonn Project",
    "D:\MPI_Brain\Data\Bonn Project",
    "C:\\Users\Albert\MyFiles\MPI_Brain\Data\Bonn Project"
    ]

for path in recognized_project_paths:
    path_exists = os.path.isdir(path)
    if path_exists:
        project_path = path
        break
else:
    raise FileNotFoundError('WARNING: No recognized project data path found in this computer')
    
    
##### MODIFIABLE PARAMETERS #####


    
# #Main data paths
#Path to the dataset
FAT_CLUSTER_PATH = project_path + "\WithAmpltd\BigFatCluster.mat"

# #Output paths
#Path to empty folders. Results store data, output stores final images.
OUTPUT_PATH = project_path + "\Output\\"
RESULTS_PATH = project_path + "\Results\\"

# #Complementary data
#Path to the place cell matlab files
PLACE_CELL_PATH_DOMBECK = project_path + "\place_cell_bool_dombeck.mat"        
PLACE_CELL_PATH_LOSONCZY = project_path + "\place_cell_bool_losonczy.mat"


# Experiment variables
MAX_POS = 1500
SESSION_NAMES = ['B1', 'B2', 'B3', 'T1', 'T2', 'Tn-1', 'Tn', 'P1', 'P2', 'Ex1','Ex2','B1','B2','T1','T2','P1','P2']


# ##### PLOTTING VARIABLES #####
PCA_CMAP = 'twilight'

 














