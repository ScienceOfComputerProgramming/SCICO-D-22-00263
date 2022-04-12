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

class AnasthesiaPKPD:

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

    def getDynamics(weight=25):
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

    def getReachSetC1C2():
        '''
        Returns the reachable set for Comparment 1 and 2, i.e., C1 and C2
        '''
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(0,6),(0,10),(0,10),(1,8),(2,2)]
        #P=[(2,4),(4,6),(4,6),(3,5),(1,3)]
        P=[(2,4),(4,6),(4,6),(3,5),(0,10)]
        initialSet=(C,V,P)

        p=1.8
        Er={
        (1,0): [1-(p/100),1+(p/100)],
        (1,1): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (2,2): [1-(p/100),1+(p/100)]
        }
        T=20

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        #exit(0)

        VisualizePKPD.vizC1C2(reachRS,reachORS)

    def getReachSetCp():
        '''
        Returns the reachable set for Comparment 1 and 2, i.e., C1 and C2
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
        initialSet=(C,V,P)

        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        T=20

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()

        #VisualizePKPD.vizCpCe(reachRS,reachORS)
        VisualizePKPD.vizCp(reachRS,reachORS)

    def getReachSetCe():
        '''
        Returns the reachable set for Comparment 1 and 2, i.e., C1 and C2
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
        initialSet=(C,V,P)

        p=5
        Er={
        (3,0): [1-(p/100),1+(p/100)],
        (3,3): [1-(p/100),1+(p/100)]
        }
        T=20

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()

        VisualizePKPD.vizCe(reachRS,reachORS)

    def onlineMonitorCp():
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P)

        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        '''p=1
        Er={
        (0,0): [1-(p/100),1+(p/100)],
        (0,4): [1-(p/100),1+(p/100)]
        }'''
        T=2000

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T)
        (l,actualBehavior)=lgr.getLog()

        #print(l[1])
        #print()
        #print(actualBehavior[l[1][0]])
        #exit(0)

        mntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSets,logs)=mntr.monitorReachSets()

        #VisualizePKPD.vizCpCe(reachRS,reachORS)
        VisualizePKPD.vizMonitorCp(reachSets,logs,T,"viz_online_monitor_cp")

    def offlineMonitorCp():
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

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P)

        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        '''p=1
        Er={
        (0,0): [1-(p/100),1+(p/100)],
        (0,4): [1-(p/100),1+(p/100)]
        }'''
        '''Er={
        (0,0): [-1,1],
        (0,4): [-1,1]
        }'''
        T=2000

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)


        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()


        mntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSets=mntr.monitorReachSets()

        #VisualizePKPD.vizCpCe(reachRS,reachORS)
        VisualizePKPD.vizMonitorCp(reachSets,logs,T,"viz_offline_monitor_cp")

    def offlineCompMonitorCp():
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        #P=[(2,4),(2,7),(3,6),(2,4),(2,10)]
        P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P)

        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        '''p=1
        Er={
        (0,0): [1-(p/100),1+(p/100)],
        (0,4): [1-(p/100),1+(p/100)]
        }'''
        '''Er={
        (0,0): [-1,1],
        (0,4): [-1,1]
        }'''
        T=2000

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)


        lgr=GenLog(A,Er,initialSet,T)
        (logs,actualBehavior)=lgr.getLog()


        mntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSets=mntr.monitorReachSets()

        print("\n--\n")

        mntr.monitorFlowStar("Anaesthesia")

        #VisualizePKPD.vizCpCe(reachRS,reachORS)
        #VisualizePKPD.vizMonitorCp(reachSets,logs,T,"viz_offline_monitor_cp")

    def compMonitorCp():
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        #P=[(2,4),(3,6),(3,6),(2,4),(2,10)]
        P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
        P_unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        P_unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafe1=(C,V,P_unsafe1)
        unsafe2=(C,V,P_unsafe2)
        unsafeList=[unsafe1,unsafe2]

        initialSet=(C,V,P)

        p=0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }
        T=2000

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        lgr=GenLog(A,Er,initialSet,T,5)
        (logs,actualBehavior)=lgr.getLog()


        onlnMntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSetsOnline,onlineLogs)=onlnMntr.monitorReachSets()

        oflnMntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSetsOffline=oflnMntr.monitorReachSets()

        VisualizePKPD.vizCompMonitorCp(reachSetsOnline,onlineLogs,reachSetsOffline,logs,T,"viz_compare_monitors_cp")










if True:
    AnasthesiaPKPD.compMonitorCp()
