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
    * Use '.' for continuous time systems (Note: MoULDyS will discretize).
    * Use '+' for discrete time systems.
    '''
    mode='+'


    ######### Optional ########
    '''
    Set discretization parameter if the mode is continuous.
    '''
    h=0.01 # This step is not needed for this example. However, having this, wouldn't have any impact in this example.


    ######### Step 3 ########
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
    * state_variable_0 >=20 AND
    * state_variable_0 <=-20
    Unsafe sets can be encoded as intervals
    '''
    unsafe1=[(-np.inf,-20),(-np.inf,np.inf),(-np.inf,np.inf)]
    unsafe2=[(20,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)]
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
        initialSet=[(0,0),(0,0),(-0.1,0.1)] # Initial Set, represented as an interval
        pr=1 # Probability of logging
        fname='toyEg'
        T=2000 # Time upto which logs are generated
        tp='interval' # Use tp='zonotope' for creating logs with zonotopes
        (log,actualBehavior)=mEngine.genLogFile(initialSet,T,fname,tp,pr)
        exit()


    ######### Step 6 ########
    '''
    Step 6: Perform offline monitoring
    * Let the log be given in file /my/location/MoULDyS/data/toyEg_20_interval (Note: Don't use .mlog extension)
    * The log type can either be interval or zonotope.

    Log file format for interval:
        * Each line: <time stamp>: <intervals>
    For example, at time step 9, say the interval is [[1,3],[3,4]], then the log file should have the following line:
        * 9: [[1,3],[3,4]]

    Log file format for zonotope:
        * Each line: <time stamp>: <center_of_zonotope>; <generator_of_zonotope>
    For example, at time step 9, say the zonotope has center=[0,0], with generator G=[[1,0],[0,1]],
    then the log file should have the following line:
        * 9: [0,0]; [[1,0],[0,1]]

    See some example log files (extension: .mlog) provided in /my/location/data/anesthesia/
    '''

    logFname='toyEg_1_interval'
    tp='interval' # Use tp='zonotope', if the log file is represented in zonotope.
    reachSets=mEngine.offlineMonitorLogFile(logFname,tp)


    ######### Step 7 ########
    '''
    Step 7: Visualize the results of monitoring
    '''
    T=2000 # Time step upto which we want to visualize
    th1=0 # State variable that is to be visualized
    vizCov=5 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
    #Note: Visualization takes time!
    mEngine.vizMonitorLogFile(reachSets,logFname,tp,T,th1,"toyEg_monitor",vizCoverage=vizCov)





toyEgOffline()
