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
from lib.VisualizeAircraft import *
from lib.ULS_Engine.StarOperations import *
from lib.GenLog import *
from lib.Monitor import *

class Aircraft:

    def createMatrix(A,B,mode,h):
        ''' Creates a single matrix based on
        . or +.
        In case of . a rough approximation is
        done'''

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

    def getDynamics(ohm=1.5):
        #ohm=1.5
        A=np.array([
        [0,0,1,0],
        [0,0,0,1],
        [0,0,0,-ohm],
        [0,0,ohm,0]
        ])
        B=np.array([])
        h=0.01
        mode='.'

        return (A,B)

    def getReachSetX1X2():
        '''
        Returns the reachable set for x1 and x2
        '''
        P=1
        Er={
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(1.1,1.11),(1.1,1.11),(20,20.1),(20,20.1)]
        initialSet=(C,V,P)
        T=400

        (dynA,dynB)=Aircraft.getDynamics(ohm=1)
        A=Aircraft.createMatrix(dynA,dynB,'.',0.01)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        print(len(reachORS))
        #exit(0)

        VisualizeAircraft.vizX1X2(reachRS,reachORS)


    def offlineMonitor():
        P=1
        Er={
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(1.1,1.11),(1.1,1.11),(20,20.1),(20,20.1)]
        initialSet=(C,V,P)
        T=2000

        P_unsafe1=[(-1000,-49),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(11,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        (dynA,dynB)=Aircraft.getDynamics(ohm=1)
        A=Aircraft.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()

        VisualizeAircraft.vizOfflineMonitorX1([],logs,T,"viz_offline_monitor_")

        mntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSets=mntr.monitorReachSets()

        VisualizeAircraft.vizOfflineMonitorX1(reachSets,logs,T,"viz_offline_monitor_")

    def onlineMonitor():
        P=0.001
        Er={
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(1.1,1.11),(1.1,1.11),(20,20.1),(20,20.1)]
        initialSet=(C,V,P)
        T=2000

        P_unsafe1=[(-1000,-49),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(11,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        (dynA,dynB)=Aircraft.getDynamics(ohm=1)
        A=Aircraft.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()

        mntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSets,logs)=mntr.monitorReachSets()

        VisualizeAircraft.vizOfflineMonitorX1(reachSets,logs,T,"viz_offline_monitor_")








if True:
    #print("Try")
    #Aircraft.getReachSetX1X2()
    #Aircraft.offlineMonitor()
    Aircraft.onlineMonitor()
