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
    Given a set of logs and a dynamics, performs monitoring.
    '''

    def __init__(self,A,Er,log,unsafe):
        self.A=A # Dynamics Matrix
        self.Er=Er # Error Matrix
        self.log=log
        self.T=len(log)
        self.unsafe=unsafe

    def monitorReachSets(self):
        '''
        Generate Reachable Sets as per logs
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
                for unsafeSet in self.unsafe:
                    int_t=0
                    for rSet in reachORS[:-1]:
                        intr=StarOp.computeIntersection(rSet,unsafeSet)
                        if intr!=-1:
                            print("\t\t>> SUBSTATUS: Computing refinement at at time step ",currentTime+int_t)
                            td=nextTime-(currentTime+int_t)
                            rsInt=Split(self.A,self.Er,intr,td)
                            (reachORSInt,reachRSInt)=rsInt.getReachableSetAllList()
                            # Check if the nextState is reachable
                            for rsi in reachORSInt[:-1]:
                                reachFlag=StarOp.checkIntersection(rsi,nextState)
                                if reachFlag==True:
                                    # nextState is reachable from unsafe set
                                    print("\t>>Safety: Unsafe at log step ",t)
                                    time_taken=time.time()-time_taken
                                    print("\t>>Time Taken: ",time_taken)
                                    return reachSets
                        int_t+=1
            else:
                if self.isSafe(reachORS[:-1])==False:
                    print("\t>>Safety: Unsafe at log step ",t)
                    time_taken=time.time()-time_taken
                    print("\t>>Time Taken: ",time_taken)
                    return reachSets
        time_taken=time.time()-time_taken
        print("\t>> Safety: Safe")
        print("\t>> Time Taken: ",time_taken)
        print(">> STATUS: Reachable sets, as per monitors, computed!")

        return reachSets

    def monitorFlowStar(self,benchmarkName):
        '''
        Apply our monitoring alorithm using Flow*
        '''
        print(">> STATUS: Using Flow* as per monitors . . .")
        time_taken=time.time()
        fs=FlowStar(benchmarkName)
        for t in range(self.T-1):
            currentTime=self.log[t][0]
            curretState=self.log[t][1]
            nextTime=self.log[t+1][0]
            nextState=self.log[t+1][1]
            timeDiff=nextTime-currentTime
            print("\t>> SUBSTATUS: Log Step: ",t,"/",self.T,"\tTime Diff: ",timeDiff)
            fg=fs.isSafe(curretState,timeDiff)
            if fg==0:
                print("\t>>Safety: Unsafe at log step ",t)
                time_taken=time.time()-time_taken
                print("\t>>Time Taken: ",time_taken)
        time_taken=time.time()-time_taken
        print("\t>> Safety: Safe")
        print("\t>> Time Taken: ",time_taken)
        print(">> STATUS: Safety using Flow*, as per monitors, computed!")


    def isSafe(self,rsList):
        '''
        Given lists of reachable sets and unsafe sets, check safety
        '''

        t=0
        for rs in rsList:
            for un in self.unsafe:
                fg=StarOp.checkIntersection(rs,un)
                if fg==True:
                    # Intersects
                    return False

        return True

    def refine(rs,t,nextLog):
        '''
        - Compute Intersection between rs and self.unsafe
        - Compute forward reachability
        - Check if the forward if the nextLog is reachable from the computed intersection
        '''



class OnlineMonitor:
    '''
    Given a dynamics, performs online monitoring.
    '''
    def __init__(self,A,Er,actualBehavior,unsafe):
        self.A=A # Dynamics Matrix
        self.Er=Er # Error Matrix
        self.actualBehavior=actualBehavior # Actual Behavior of the system
        self.T=len(actualBehavior)
        self.unsafe=unsafe

    def monitorReachSets(self):
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
                    print("\t>> Safety: Truly unsafe at time ",t+1)
                    time_taken=time.time()-time_taken
                    print("\t>>Time Taken: ",time_taken)
                    return (reachSets,logs)
            reachSets.append(copy.copy(reachSet))
        time_taken=time.time()-time_taken
        print("\t>> Safety: Safe")
        print("\t>> Number of Logs: ",len(logs))
        print("\t>> Time Taken: ",time_taken)
        print(">> STATUS: Performed online monitoring using reachable sets!")
        return (reachSets,logs)

    def monitorFlowStar(self,benchmarkName):
        # Triggering a log at the first step
        reachSet=copy.copy(self.triggerLog(0))
        reachSets=[copy.copy(reachSet)]
        logs=[[0,copy.copy(reachSet)]]
        fs=FlowStar(benchmarkName)
        print(">> STATUS: Performing online monitoring using Flow* . . .")
        time_taken=time.time()
        for t in range(self.T-1):
            # Compute one step reachable set
            rs=Split(self.A,self.Er,reachSet,1)
            (reachORS,reachRS)=rs.getReachableSetAllList()
            reachSet=copy.copy(reachORS[1])
            # Trigger a log if unsafe at t+1
            if self.isSafe([reachSet])==False:
                print("\t>> SUBSTATUS: Triggering a log at time ",t+1)
                # Trigger a log
                reachSet=copy.copy(self.triggerLog(t+1))
                logs.append([t+1,reachSet])
                if self.isSafe([reachSet])==False:
                    print("\t>> Safety: Truly unsafe at time ",t+1)
                    time_taken=time.time()-time_taken
                    print("\t>>Time Taken: ",time_taken)
                    return (reachSets,logs)
            reachSets.append(copy.copy(reachSet))
        time_taken=time.time()-time_taken
        print("\t>> Safety: Safe")
        print("\t>> Number of Logs: ",len(logs))
        print("\t>> Time Taken: ",time_taken)
        print(">> STATUS: Performed online monitoring using Flow*!")
        return (reachSets,logs)

    def triggerLog(self,t):
        return self.actualBehavior[t]

    def isSafe(self,rsList):
        '''
        Given lists of reachable sets and unsafe sets, check safety
        '''

        for rs in rsList:
            for un in self.unsafe:
                fg=StarOp.checkIntersection(rs,un)
                if fg==True:
                    # Intersects
                    return False

        return True
