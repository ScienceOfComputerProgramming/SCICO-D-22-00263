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
from lib.Visualize import *
from lib.ULS_Engine.StarOperations import *
from lib.GenLog import *
from lib.Monitor import *

'''
This is to just test installation.

Run this script: If it runs without error, the environment is ready!
'''


class TestInstall:

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

    def tester():
        print(f"{bcolors.OKGREEN}>> STATUS: Started testing environment . . .{bcolors.ENDC}")
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        P=[(15,15.01),(3,3.03),(14.9,15),(1,1)]
        P_unsafe=[(-1000,1000),(-1000,0.5),(-1000,1000),(-1000,1000),(-1000,1000)]
        unsafe=(C,V,P_unsafe)
        unsafeList=[unsafe]
        initialSet=(C,V,P)
        Er={
        (0,3): [-0.6,2.46],
        (2,3): [-0.9,0.6]
        }
        T=10
        (dynA,dynB)=TestInstall.getDynamics()
        A=TestInstall.createMatrix(dynA,dynB,'.',0.01)
        lgr=GenLog(A,Er,initialSet,T,14)
        (logs,actualBehavior)=lgr.getLog()
        onlnMntr=OnlineMonitor(A,Er,actualBehavior,unsafeList)
        (reachSetsOnline,onlineLogs)=onlnMntr.monitorReachSets()
        oflnMntr=OfflineMonitor(A,Er,logs,unsafeList)
        reachSetsOffline=oflnMntr.monitorReachSets()
        Visualize.vizTest(reachSetsOnline,onlineLogs,T,"viz_test")
        print(f"{bcolors.OKGREEN}>> STATUS: Environment testing completed!{bcolors.ENDC}")
        print(f"\n\n{bcolors.OKCYAN}{bcolors.BOLD}=======================")
        print("Environment is Ready!")
        print(f"======================={bcolors.ENDC}\n\n")


if True:
    TestInstall.tester()
