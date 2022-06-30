'''
Implements online and offline monitoring algorithms as per [1]

[1] Bineet Ghosh and Étienne André.
    "Offline and online monitoring of scattered uncertain logs using uncertain linear dynamical systems."
    In: FORTE 2022.
'''

import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from Parameters import *
import random
from lib.ULS_Engine.SplitMet import *
from lib.ULS_Engine.StarOperations import *
#from lib.FlowStar import *
import copy
import time

class OfflineMonitor:
    '''
    Given a set of logs and a dynamics, performs monitoring, as per Section 4.1 of the paper.
    In other words, implements Algorithm 1 of the paper.
    '''

    def __init__(self,A,Er,log,unsafe):
        self.A=A # Dynamics Matrix
        self.Er=Er # Error Matrix
        self.log=log # The given log as input
        self.T=len(log) # The maximum time horizon
        self.unsafe=unsafe # The given unsafe specification

    def monitorReachSets(self):
        '''
        Implements the offline monitoring, Algorithm 1, as given in [1].

        Output:
            1. Prints the safety result: Potentially-unsafe / Safe.
            2. Returns the computed reachable  sets.
        '''

        reachSets=[]

        print(">> STATUS: Computing reachable sets as per monitors . . .")
        time_taken=time.time()
        for t in range(self.T-1):
            currentTime=self.log[t][0]
            curretState=self.log[t][1]
            nextTime=self.log[t+1][0]
            nextState=self.log[t+1][1]
            timeDiff=nextTime-currentTime
            print("\t>> SUBSTATUS: Log Step: ",t,"/",self.T,"\tTime Diff: ",timeDiff)
            rs=Split(self.A,self.Er,curretState,timeDiff)
            (reachORS,reachRS)=rs.getReachableSetAllList()
            reachSets.extend(reachORS[:-1])
            if REFINE==True:
                # The refinement module starts.
                for unsafeSet in self.unsafe:
                    int_t=0
                    for rSet in reachORS[:-1]:
                        intr=StarOp.computeIntersection(rSet,unsafeSet)
                        if intr!=-1:
                            print("\t\t>> SUBSTATUS: Computing refinement at at time step ",currentTime+int_t)
                            td=nextTime-(currentTime+int_t)
                            rsInt=Split(self.A,self.Er,intr,td)
                            (reachORSInt,reachRSInt)=rsInt.getReachableSetAllList() # Computes reach set of uncertain linear systems
                            # Check if the nextState is reachable
                            for rsi in reachORSInt[:-1]:
                                reachFlag=StarOp.checkIntersection(rsi,nextState)
                                if reachFlag==True:
                                    # nextState is reachable from unsafe set
                                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Safety:{bcolors.ENDC} {bcolors.BOLD}{bcolors.FAIL}Unsafe at log step ",t,f"{bcolors.ENDC}")
                                    time_taken=time.time()-time_taken
                                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Time Taken: ",time_taken,f"{bcolors.ENDC}")
                                    return reachSets
                        int_t+=1
            else:
                if self.isSafe(reachORS[:-1])==False:
                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Safety:{bcolors.ENDC} {bcolors.BOLD}{bcolors.FAIL}Unsafe at log step ",t,f"{bcolors.ENDC}")
                    time_taken=time.time()-time_taken
                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Time Taken: ",time_taken,f"{bcolors.ENDC}")
                    return reachSets
        time_taken=time.time()-time_taken
        print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Safety:{bcolors.ENDC} {bcolors.BOLD}{bcolors.OKGREEN}Safe{bcolors.ENDC}")
        print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Time Taken: ",time_taken,f"{bcolors.ENDC}")
        print(">> STATUS: Reachable sets, as per monitors, computed!")

        return reachSets

    def isSafe(self,rsList):
        '''
        Given lists of reachable sets and unsafe sets, check safety.

        Returns true iff none of elements in `rsList` intersects with the unsafe set.
        '''

        t=0
        for rs in rsList:
            for un in self.unsafe:
                fg=StarOp.checkIntersection(rs,un)
                if fg==True:
                    # Intersects
                    return False

        return True




class OnlineMonitor:
    '''
    Given a dynamics, performs online monitoring, as per Section 4.2 of [1].
    In other words, implements Algorithm 2 of [1].
    '''
    def __init__(self,A,Er,actualBehavior,unsafe):
        self.A=A # Dynamics Matrix
        self.Er=Er # Error Matrix
        self.actualBehavior=actualBehavior # Actual Behavior of the system
        self.T=len(actualBehavior) # The maximum time horizon
        self.unsafe=unsafe # The given unsafe specification

    def monitorReachSets(self):
        '''
        Given a dynamics, performs online monitoring, as per Section 4.2 of [1].
        In other words, implements Algorithm 2 of [1].

        Output:
            1. Prints the safety result: Potentially-unsafe / Safe.
            2. Returns the computed reachable  sets and the logs.
        '''
        # Triggering a log at the first step
        reachSet=copy.copy(self.triggerLog(0))
        time_since_last_log=time.time()
        reachSets=[copy.copy(reachSet)]
        logs=[[0,copy.copy(reachSet)]]
        print(">> STATUS: Performing online monitoring using reachable sets . . .")
        time_taken=time.time()
        for t in range(self.T-1):
            # Compute one step reachable set
            rs=Split(self.A,self.Er,reachSet,1)
            (reachORS,reachRS)=rs.getReachableSetAllList()
            reachSet=copy.copy(reachORS[1])
            # Trigger a log if unsafe at t+1
            if self.isSafe([reachSet])==False:
                print("\t>> SUBSTATUS: Triggering a log at time ",t+1,"\t Time Taken (Diff): ",time.time()-time_since_last_log)
                # Trigger a log
                reachSet=copy.copy(self.triggerLog(t+1))
                time_since_last_log=time.time()
                logs.append([t+1,reachSet])
                if self.isSafe([reachSet])==False:
                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Safety:{bcolors.ENDC} {bcolors.BOLD}{bcolors.FAIL}Truly unsafe at time ",t+1,f"{bcolors.ENDC}")
                    time_taken=time.time()-time_taken
                    print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Time Taken: ",time_taken,f"{bcolors.ENDC}")
                    return (reachSets,logs)
            reachSets.append(copy.copy(reachSet))
        time_taken=time.time()-time_taken
        print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Safety:{bcolors.ENDC} {bcolors.BOLD}{bcolors.OKGREEN}Safe{bcolors.ENDC}")
        print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Number of Logs: ",len(logs),f"{bcolors.ENDC}")
        print(f"\t{bcolors.OKCYAN}{bcolors.BOLD}>> Time Taken: ",time_taken,f"{bcolors.ENDC}")
        print(">> STATUS: Performed online monitoring using reachable sets!")
        return (reachSets,logs)

    def triggerLog(self,t):
        '''
        Given a time stamp `t`, returns the actual behavior at that time step.
        '''
        return self.actualBehavior[t]

    def isSafe(self,rsList):
        '''
        Given lists of reachable sets and unsafe sets, check safety.

        Returns true iff none of elements in `rsList` intersects with the unsafe set.
        '''

        for rs in rsList:
            for un in self.unsafe:
                fg=StarOp.checkIntersection(rs,un)
                if fg==True:
                    # Intersects
                    return False

        return True
