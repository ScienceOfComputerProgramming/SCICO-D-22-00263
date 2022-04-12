import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from Parameters import *
import random
from lib.ULS_Engine.ComputeU import *
import copy

class GenLog:
    '''
    Given a linear uncertain system and an initial set, generate logs for monitoring
    '''
    def __init__(self,A,Er,initialSet,T,pr=PROBABILITY_LOG):
        self.A=A # The nominal dynamics matrix
        self.Er=Er # Uncertainties
        self.initialSet=initialSet # Initial Set
        self.T=T # The maximum time horizon
        self.pr=pr # The probability of logging at a given step is p/100

    def getAccLog(self):
        '''
        Generate a log:
            - Take a single point from initial set.
            - Assume a fixed valuation of the uncertainties.
            - At each step, log values with probability p.
        '''
        log=[] # log: a list of reachable points.

        '''# Sample a random point from the initial set
        C_p=self.initialSet[0]
        V_p=self.initialSet[1]
        P_p=[]
        for p in self.initialSet[2]:
            pt=random.uniform(p[0],p[1])
            P_p.append((pt,pt))
        initPoint=(C_p,V_p,P_p)'''
        initPoint=self.initialSet
        log.append([0,initPoint])
        actualBehavior=[initPoint]

        # Reachable Set of A_point starting from initPoint
        PROB=100
        rs=copy.copy(initPoint)
        for t in range(self.T-1):
            # Sample a random uncertaity
            A_point=copy.copy(self.A)
            for key in self.Er:
                A_point[key[0]][key[1]]=self.A[key[0]][key[1]]*random.uniform(self.Er[key][0],self.Er[key][1])
            rs=CompU.prodMatStars(A_point,rs)
            # To log, or to not
            rdBit=random.randint(1,PROB)
            if rdBit<=self.pr:
                # Log
                log.append([t+1,rs])
            actualBehavior.append(rs)

        return (log,actualBehavior)

    def getLog(self):
        '''
        Get a log with uncertaity in it:
            - Generate an accurate log.
            - Add uncertaity to the reachable set.
            - Add uncertaity to time. (ToDo)
        '''
        print(">>STATUS: Generating a log  . . .")
        time_taken=time.time()
        accLog,actualBehavior=self.getAccLog()
        log=[]
        actualBehaviorUn=[]

        C=[0]*self.A.shape[0]
        V=np.identity(self.A.shape[0])
        P=[(-0.4,0.4)]*self.A.shape[0]
        unitBox=(C,V,P)

        fg=0
        for lg in accLog:
            if fg!=0:
                rsLg=CompU.addStars(lg[1],unitBox) # Additing uncertaity to reachable set
                log.append([lg[0],rsLg])
            else:
                log.append([lg[0],lg[1]])
            fg=1

        fg=0
        for lg in actualBehavior:
            if fg!=0:
                rsLg=CompU.addStars(lg,unitBox) # Additing uncertaity to reachable set
                actualBehaviorUn.append(rsLg)
            else:
                actualBehaviorUn.append(lg)
            fg=1

        print("\t>>STATUS: Time Taken: ",time.time()-time_taken)
        print(">>STATUS: Log generated!")

        return (log,actualBehaviorUn)


    def getLogs(self,N=5):
        logs=[]
        print(">>STATUS: Generating logs  . . .")
        time_taken=time.time()
        for i in range(N):
            logs.append(self.getLog())
        print("\t>>STATUS: Time Taken: ",time.time()-time_taken)
        print(">>STATUS: Logs generated!")

        return logs
