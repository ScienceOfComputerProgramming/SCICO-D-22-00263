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


def toyEgOnline():
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



    ######### Step 6 ########
    '''
    Step 6: Perform online monitoring
    * Let the actual behavior be given in file /my/location/MoULDyS/data/toyEg_1_interval (Note: Don't use .mlog extension)
    * The actual behavior type can either be interval or zonotope.

    See some example log files (extension: .mbeh) provided in /my/location/data/anesthesia/
    '''

    logFname='toyEg_1_interval'
    tp='interval' # Use tp='zonotope', if the log file is represented in zonotope.
    (reachSets,logs)=mEngine.onlineMonitorBehFile(logFname,tp)


    ######### Step 7 ########
    '''
    Step 7: Visualize the results of monitoring
    '''
    T=2000 # Time step upto which we want to visualize
    th1=0 # State variable that is to be visualized
    vizCov=5 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
    #Note: Visualization takes time!
    mEngine.vizMonitor(reachSets,logs,tp,T,th1,"toyEg_monitor",vizCoverage=vizCov)





toyEgOnline()
