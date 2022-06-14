import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from Parameters import *

import numpy as np
import math

from lib.ULS_Engine.SplitMet import *
from lib.ULS_Engine.OrderUncertainties import *
from lib.VisualizePKPD import *
from lib.ULS_Engine.StarOperations import *
from lib.GenLog import *
from lib.Monitor import *

'''
This benchmark has been taken from [3].
Illustrates online and offline monitoring algorithms, as per [1], on the example from [3].


[1] Bineet Ghosh and Étienne André.
    "Offline and online monitoring of scattered uncertain logs using uncertain linear dynamical systems."
    In: FORTE 2022.
[3] Victor Gan, Guy Albert Dumont, and Ian Mitchell. “Benchmark Problem: A PK/PD Model and Safety Constraints for Anesthesia Delivery”.
    In: Proceedings of the 1st and 2nd International Workshops on Applied veRification for Continuous and Hybrid Systems.
'''

class AnasthesiaPKPD:

    @staticmethod
    def createMatrix(A,B,mode,h):
        ''' Creates a single matrix based on
        `.` or `+`. Here, `.` represents contionous systems, `+` represents discrete system.

        In case of `.` a rough approximation is performed.

        Return: This function augments the matrices A and B to form a single matrix.
        '''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n),dtype=np.float)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def getDynamics(weight=25):
        '''
        The dynamics of the benchmark as per [3].

        Returns: The dynamics matrices A and B.
        '''
        #With weight of 25 Kg
        v1=458.4*weight
        k10=0.1527*pow(weight,-0.3)
        k12=0.114
        k13=0.0419
        k21=0.055
        k31=0.0033
        kd=40 #with Td=20s
        A=np.array([
        [-(k10+k12+k13),k12,k13,0],
        [k21,-k21,0,0],
        [k31,0,-k31,0],
        [kd,0,0,-kd]
        ])
        B=np.array([
        [1/v1],
        [0],
        [0],
        [0]
        ])

        return (A,B)

    def onlineMonitorCp():
        '''
        Performs offline monitoring of c_p (the plasma concentration level).
        '''
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        #P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        '''
        The behavior is marked as unsafe if: c_p<=1 OR c_p>=6. This encoded by `unsafe`.
        '''
        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P)

        '''
        Uncertainties in the given model of PKPD
        '''
        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }

        T=2000 # Monitoring will be performed upto time step `T`.

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01) # Augment the matrices A and B to form a single matrix.


        '''
        Obtain a simulated actual behavior. In case the use already has access to
        the actual behavior, this step is not needed---they can instead use the
        available/desired behavior. This is needed for the logger to get access
        to the actual behavior when needed
        '''
        lgr=GenLog(A,Er,initialSet,T)
        (l,actualBehavior)=lgr.getLog()

        mntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSets,logs)=mntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])

        # Note: This step takes quite some time.
        VisualizePKPD.vizMonitorCp(reachSets,logs,T,"viz_online_monitor_cp") # Visualize the logs and the reachable sets.

    def offlineMonitorCp():
        '''
        Performs online monitoring of c_p (the plasma concentration level).
        '''
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        #P=[(2,4),(2,7),(3,6),(2,4),(2,10)]
        #P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        '''
        The behavior is marked as unsafe if: c_p<=1 OR c_p>=6. This encoded by `unsafe`.
        '''
        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P) # This is the initial set.

        '''
        Uncertainties in the given model of PKPD
        '''
        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }

        T=2000 # Monitoring will be performed upto time step `T`.

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01) # Augment the matrices A and B to form a single matrix.

        '''
        Obtain a simulated log. In case the use already has access to a log,
        this step is not needed---they can instead use the available/desired log
        that is to be monitored.
        '''
        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()


        mntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSets=mntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])

        # Note: This step takes quite some time.
        VisualizePKPD.vizMonitorCp(reachSets,logs,T,"viz_offline_monitor_cp") # Visualize the logs and the reachable sets.

    def compMonitorCp():
        '''
        Similar to the above two functions, this API compares the online and offline monitoring.
        '''
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        #P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        '''
        The behavior is marked as unsafe if: c_p<=1 OR c_p>=6. This encoded by `unsafe`.
        '''
        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P) # This is the initial set.

        '''
        Uncertainties in the given model of PKPD
        '''
        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        T=2000 # Monitoring will be performed upto time step `T`.

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01) # Augment the matrices A and B to form a single matrix.

        '''
        Obtain a simulated actual behavior. In case the use already has access to
        the actual behavior, this step is not needed---they can instead use the
        available/desired behavior. This is needed for the logger to get access
        to the actual behavior when needed. And the `logs` is needed by the offline
        monitoring algorithm as an input.
        '''
        lgr=GenLog(A,Er,initialSet,T,5)
        (logs,actualBehavior)=lgr.getLog()

        onlnMntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSetsOnline,onlineLogs)=onlnMntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])

        oflnMntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSetsOffline=oflnMntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])

        # Note: This step takes quite some time.
        # This visualizes the results from both online and offline monitoring.
        VisualizePKPD.vizCompMonitorCp(reachSetsOnline,onlineLogs,reachSetsOffline,logs,T,"viz_compare_monitors_cp")










if True:
    AnasthesiaPKPD.offlineMonitorCp()
