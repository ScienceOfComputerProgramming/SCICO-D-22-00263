import pickle

import sys,os

import numpy as np
import math

from Parameters import *
from SplitMet import *
from OrderUncertainties import *
from VisualizePKPD import *
from StarOperations import *

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

    def getCellOrder():
        (A,B)=AnasthesiaPKPD.getDynamics()
        A_club=AnasthesiaPKPD.createMatrix(A,B,'.',0.01)
        print(">> STATUS: Applying Cell Ordering Algorithm . . .")
        time_taken=time.time()
        projMat=np.zeros(A_club.shape)
        projMat[1][1]=1
        projMat[0][0]=0
        projMat[2][2]=0
        A_club_proj=np.matmul(projMat,A_club)
        cellOrder=OrdUnc(A_club_proj).getOrder()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Applying Cell Ordering Algorithm . . .")
        print(cellOrder)

    def getTopBotErrors():
        '''
        Based on the ouput we got from `self.getCellOrder()`,
        some perturbation of the matrix were picked.
        '''
        p=2

        Er_top={
        (1,0): [1-(p/100),1+(p/100)],
        (1,1): [1-(p/100),1+(p/100)]
        }

        Er_bot={
        (2,0): [1-(p/100),1+(p/100)],
        (2,2): [1-(p/100),1+(p/100)]
        }

        return Er_top,Er_bot

    def getReachSetsTopBot():
        '''
        Based on Errors obtained from cell-ordering algorithm, we
        compare the top and the bottom cells
        '''
        (Er_top,Er_bot)=AnasthesiaPKPD.getTopBotErrors()
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(2,4),(4,6),(4,6),(3,5),(0,10)]
        initialSet=(C,V,P)
        T=20

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        rsTop=Split(A,Er_top,initialSet,T)
        (reachORSTop,reachRSTop)=rsTop.getReachableSetAllList()

        rsBot=Split(A,Er_bot,initialSet,T)
        (reachORSBot,reachRSBot)=rsBot.getReachableSetAllList()

        VisualizePKPD.vizComp(reachORSTop,reachORSBot)

    def getRobustnessMetric():
        '''
        Returns the Robustness Metric
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
        C_u=[0,0,0,0,0]
        V_u=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P_u=[(6,60),(-30000000000000,60000000000000),(-30000000000000,60000000000000),(-30000000000000,60000000000000),(-30000000000000,60000000000000)]
        unsafe1=(C_u,V_u,P_u)
        P_u_2=[(-60,1.2),(-30000000000000,60000000000000),(-30000000000000,60000000000000),(-30000000000000,60000000000000),(-30000000000000,60000000000000)]
        unsafe2=(C_u,V_u,P_u_2)

        p=0.1
        T=20

        (dynA,dynB)=AnasthesiaPKPD.getDynamics()
        A=AnasthesiaPKPD.createMatrix(dynA,dynB,'.',0.01)

        step=0.1
        rsListOld=[]
        rsLisUnpOld=[]
        time_taken=time.time()
        while p<=0.8:
            Er={
            (0,0): [1-(p/2000),1+(p/2000)],
            (0,4): [100/(100+p),100/(100-p)]
            }
            rs=Split(A,Er,initialSet,T)
            (reachORS,reachRS)=rs.getReachableSetAllList()
            safe=True
            for rs in reachORS:
                int1=StarOp.checkIntersection(rs,unsafe1)
                int2=StarOp.checkIntersection(rs,unsafe2)
                if int1==True or int2==True:
                    safe=False
                    p=p-step
                    break;
            if safe==False:
                break
            p=p+step
            rsListOld=copy.copy(reachORS)
            rsListUnpOld=copy.copy(reachRS)
        time_taken=time.time()-time_taken
        print("\tRobustness Metric: ",p)
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Robustness Metric Computed")
        print(">> STATUS: Visualizing Safe Reachable Sets  . . .")
        VisualizePKPD.vizCp(rsListUnpOld,rsListOld,fname="viz_Cp_robMet")






if True:
    AnasthesiaPKPD.getRobustnessMetric()
