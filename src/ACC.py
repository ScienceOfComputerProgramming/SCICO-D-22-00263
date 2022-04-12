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
from lib.VisualizeACC import *
from lib.ULS_Engine.StarOperations import *
from lib.GenLog import *
from lib.Monitor import *

class ACC:

    @staticmethod
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

    def getDynamics():
        m=1370
        g=9.8
        f1=1.2567
        F=-1.2156
        aL=0.5
        A=np.array([
        [f1/m,0,0,F],
        [-1,0,1,0],
        [0,0,0,aL],
        [0,0,0,0]
        ])
        B=np.array([])
        h=0.01
        mode='.'

        return (A,B)


    def offlineMonitorH():
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        #P=[(15,15.01),(3,3.03),(14.9,15),(1,1)]
        P=[(15,15.1),(3,3.5),(14.9,15.1),(1,1)]
        #P=[(15,16.5),(4,6.5),(14,15.5),(1,1)]


        P_unsafe=[(-1000,1000),(-1000,0.5),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe=(C,V,P_unsafe)
        unsafeList=[unsafe]

        initialSet=(C,V,P)

        Er={
        (0,3): [-0.6,2.46],
        (2,3): [-0.9,0.2]
        }
        T=2000

        (dynA,dynB)=ACC.getDynamics()
        A=ACC.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()


        mntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSets=mntr.monitorReachSets()

        VisualizeACC.vizMonitorH(reachSets,logs,T,"viz_offline_monitor_h")

    def onlineMonitorH():
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        #P=[(16,22),(20,25),(16,22),(1,1)]
        #P=[(15,24),(18,27),(15,24),(1,1)]
        #P=[(20,22),(20,23),(20,22),(1,1)]
        P=[(15,15.01),(3,3.03),(14.9,15),(1,1)]

        P_unsafe=[(-1000,1000),(-1000,0.5),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe=(C,V,P_unsafe)
        unsafeList=[unsafe]

        initialSet=(C,V,P)

        Er={
        (0,3): [-0.6,2.46],
        (2,3): [-0.9,0.6]
        }
        T=2000

        (dynA,dynB)=ACC.getDynamics()
        A=ACC.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T)
        (l,actualBehavior)=lgr.getLog()


        mntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSets,logs)=mntr.monitorReachSets()

        #VisualizePKPD.vizCpCe(reachRS,reachORS)
        VisualizeACC.vizMonitorH(reachSets,logs,T,"viz_online_monitor_h")


    def compMonitorH():
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        #P=[(16,22),(20,25),(16,22),(1,1)]
        #P=[(15,24),(18,27),(15,24),(1,1)]
        #P=[(20,22),(20,23),(20,22),(1,1)]
        P=[(15,15.01),(3,3.03),(14.9,15),(1,1)]

        P_unsafe=[(-1000,1000),(-1000,0.5),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe=(C,V,P_unsafe)
        unsafeList=[unsafe]

        initialSet=(C,V,P)

        Er={
        (0,3): [-0.6,2.46],
        (2,3): [-0.9,0.6]
        }
        T=2000

        (dynA,dynB)=ACC.getDynamics()
        A=ACC.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T,14)
        (logs,actualBehavior)=lgr.getLog()


        onlnMntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSetsOnline,onlineLogs)=onlnMntr.monitorReachSets()

        oflnMntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSetsOffline=oflnMntr.monitorReachSets()

        VisualizeACC.vizCompMonitorH(reachSetsOnline,onlineLogs,reachSetsOffline,logs,T,"viz_compare_monitors_h")


if False:
    ACC.onlineMonitorH()
