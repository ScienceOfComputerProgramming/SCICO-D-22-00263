'''
This is a hello world file for MoULDyS to perform offline monitoring.

We use a toy dynamics as follows:
A=[
[1, c],
[0, 1]
],
c \in 0.1 +/- 1%,
B=[
[d],
[0.01]
]
d \in 0.01 +/- 1%,
x[t+1]=Ax[t]+BU,

where, x[0]=[[0,0],[0,0]] and U \in [[-0.1,0.1]].
'''

'''
Step 0.

Set the project path and import MoULDyS engine.
'''

import os,sys
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT) # Set the project path

from lib.MoULDySEngine import * # Importing all the functionalities of MoULDyS.


def toyEgOffline():
    ######### Step 1 ########
    '''
    Step 1: Define the dynamics matrices
    '''
    A=np.array(
    [
    [1, 0.1],
    [0, 1]
    ]
    )
    B=np.array(
    [
    [0.01],
    [0.01]
    ]
    )


    ######### Step 2 ########
    '''
    Step 2: Define the mode of the dynamics.
    * Use '.' for contionous time systems (Note: MoULDyS will discretize).
    * Use '+' for discrete time systems.
    '''
    mode='+'


    ######### Optional ########
    '''
    Set discretization parameter if the mode is contionous.
    '''
    h=0.01 # This step is not needed for this example. However, having this, wouldn't have any impact in this example.


    ######### Step 2 ########
    '''
    Step 3: Encode the uncertainties in the dynamics (i.e, c and d)
    '''
    Er={
    (0,1): [0.99,1.01], # Encoding c
    (0,2): [0.99,1.01], # Encoding d
    }


    ######### Step 4 ########
    '''
    Step 4: Encode the unsafe set.
    Let the unsafe behavior be as follows:
    * state_variable_0 >=200 AND
    * state_variable_0 <=-200
    Unsafe sets can be encoded as intervals
    '''
    unsafe1=[(-np.inf,-200),(-np.inf,-np.inf),(-np.inf,-np.inf)]
    unsafe2=[(200,np.inf),(np.inf,np.inf),(np.inf,np.inf)]
    unsafeList=[unsafe1,unsafe2]


    ######### Step 5 ########
    '''
    Step 5: Instantiate the MoULDyS engine.
    '''
    mEngine=MoULDyS(A,B,Er,mode,unsafeList,h) # Note: In this example, h is optional.


    ######### Optinal ########
    '''
    (This step can be skipped) Set the flag to true if you want to generate random logs.
    '''
    if False:
        initialSet=[(0,0),(0,0),(-0.1,0.1)]
        pr=40 # Probability of logging
        fname='toyEg'
        T=2000 # Time upto which logs are generated
        tp='zonotope'
        (log,actualBehavior)=mEngine.genLogFile(initialSet,T,fname,tp,pr)
        exit()


    ######### Step 6 ########
    '''
    Step 6: Perform offline monitoring
    * Let the log be given in file /my/location/MoULDyS/data/toyEg_20_interval (Note: Don't use .mlog extension)
    * The log type can either be interval or zonotope.

    Log file format for interval:
        * Each line: <time stamp>: <intervals>

    Log file format for zonotope:
        * Each line: <time stamp>: <center_of_zonotope>; <generator_of_zonotope>

    See some example log files (extension: .mlog) provided in /my/location/data/anesthesia/
    '''

    logFname='toyEg_40_interval'
    tp='interval' # Use tp='zonotope', if the log file is represented in zonotope.
    reachSets=mEngine.offlineMonitorLogFile(logFname,tp)


    ######### Step 7 ########
    '''
    Step 7: Visualize the results of monitoring
    '''
    # ToDo





toyEgOffline()
