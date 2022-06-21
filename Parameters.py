'''
All the parameters required for the tool can be set from this file.
'''


import os,sys

'''
Please add the following line in ~/.bashrc
export MNTR_ROOT_DIR = <YOUR PROJECT ROOT>
'''

PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

LIB_PATH=PROJECT_ROOT+'/'+'lib/' # Path to all the library functions, including the main algorithms.
SRC_PATH=PROJECT_ROOT+'/'+'src/' # Path to all the case study and example implementations.
OUTPUT_PATH=PROJECT_ROOT+'/'+'output/' # Path where all the plots are saved.
PICKLE_PATH=PROJECT_ROOT+'/'+'pickles/' # Path where all the pickle files are stored.
DATA_PATH=PROJECT_ROOT+'/'+'data/' # Path where all the data files (inputs) are stored.

'''
Define colors for colored inputs.
'''
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


PICKLE_FLAG=True # Set this flag to true, if the monitoring results are to be pickled.
REFINE=True # Set this flag to true, if you want to execute the refinement module in offline monitoring (Algorithm 2).
# Note: We strongly recommend REFINE=True


'''
Visualizing the reachable sets take a lot of time. We do not recommend visulaizing all reachable sets.
Visualizing some reachable sets, periodically, has been observed to sufficient.
'''
VIZ_PER_COVERAGE=20 # Percentage of reachable sets to be visualized.
PROBABILITY_LOG=1 # Percentage probability of logging at each time step.
