'''
A Wrapper around `Monitor.py`.

Eases the process to perform Offline and Online Monitoring
'''
import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from Parameters import *
from lib.Monitor import *
from lib.InpParser import *
from lib.Visualize import *

class MoULDyS:

    def __init__(self,A,B,Er,mode,unsafeList,h=0.01):
        self.A=MoULDyS.createMatrix(A,B,mode,h)
        self.Er=Er # Error Matrix
        self.unsafeList=self.castUnsafety(unsafeList) # The given unsafe specification

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
            I=np.identity(n1,dtype=np.float)
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

    def castUnsafety(self,unsafeList):
        mUnsafeList=[]
        for un in unsafeList:
            C=[0]*self.A.shape[0]
            V=np.identity(self.A.shape[0])
            mUnsafeList.append((C,V,un))
        return mUnsafeList

    def castIS(self,rs):
        mRS=[]
        C=[0]*self.A.shape[0]
        V=np.identity(self.A.shape[0])
        mRS=(C,V,rs)
        return mRS


    def offlineMonitor(self,logs):
        mntr=OfflineMonitor(self.A,self.Er,logs,self.unsafeList)
        reachSets=mntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])
        return reachSets

    def offlineMonitorLogFile(self,logFname,tp='interval'):
        '''
        tp: Representation type of the log
        '''
        prsr=InpParse(logFname,tp)
        logs=prsr.getLog()
        mntr=OfflineMonitor(self.A,self.Er,logs,self.unsafeList)
        reachSets=mntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])
        return reachSets

    def onlineMonitor(self,actualBehavior):
        mntr=OnlineMonitor(self.A,self.Er,actualBehavior,self.unsafeList)
        (reachSets,logs)=mntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])
        return (reachSets,logs)

    def onlineMonitorBehFile(self,logFname,tp='interval'):
        # tp: Representation type of the log
        prsr=InpParse(logFname,tp)
        beh=prsr.getBehavior()
        mntr=OnlineMonitor(self.A,self.Er,beh,self.unsafeList)
        (reachSets,logs)=mntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])
        return (reachSets,logs)

    def compMonitor(self,logs,actualBehavior):
        # actualBehavior: Actual behavior of the system
        mntr=OfflineMonitor(self.A,self.Er,logs,self.unsafeList)
        reachSetsOffline=mntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])

        mntr=OnlineMonitor(self.A,self.Er,actualBehavior,self.unsafeList)
        (reachSetsOnline,logsOnline)=mntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])

        return (reachSetsOffline,reachSetsOnline,logsOnline)

    def compMonitorFile(self,logFname,tp='interval'):
        '''
        tp: Representation type of the log
        '''
        prsr=InpParse(logFname,tp)
        logs=prsr.getLog()
        beh=prsr.getBehavior()

        mntr=OfflineMonitor(self.A,self.Er,logs,self.unsafeList)
        reachSetsOffline=mntr.monitorReachSets() # Perform offline monitoring (Algorithm 1 of [1])

        mntr=OnlineMonitor(self.A,self.Er,actualBehavior,self.unsafeList)
        (reachSetsOnline,logsOnline)=mntr.monitorReachSets() # Perform online monitoring (Algorithm 2 of [1])

        return (reachSetsOffline,reachSetsOnline,logsOnline)

    def genLogFile(self,initialSetInt,T,fname,tp,pr=PROBABILITY_LOG):
        '''
        initialSetInt: Initial set
        T: Max time step
        tpRep: Representation type of the log
        pr: Probability of logging
        '''
        initialSet=self.castIS(initialSetInt)
        lgr=GenLog(self.A,self.Er,initialSet,T,pr)
        (l,actualBehavior)=lgr.getLogFile(fname+"_"+str(pr)+"_"+tp,tp)
        #print(actualBehavior[0])
        return (l,actualBehavior)

    def genLog(self,initialSetInt,T,pr=PROBABILITY_LOG):
        '''
        initialSetInt: Initial set
        T: Max time step
        tpRep: Representation type of the log
        pr: Probability of logging
        '''
        initialSet=self.castIS(initialSetInt)
        lgr=GenLog(self.A,self.Er,initialSet,T,pr)
        (l,actualBehavior)=lgr.getLog()
        #print(actualBehavior[0])
        return (l,actualBehavior)

    def vizMonitor(self,ORS_List,logs,tp,T,th1,fname="viz_test",vizCoverage=VIZ_PER_COVERAGE):
        '''
        ORS_List: List of computed reachable sets
        logs: Input logs to monitor
        tp: Representation type of the log
        th1: The state variable to be visualized
        T: Max time step
        '''
        Visualize.vizMonitor(ORS_List,logs,T,th1,self.unsafeList,fname,vizCoverage)

    def vizMonitorLogFile(self,ORS_List,logFname,tp,T,th1,fname="viz_test",vizCoverage=VIZ_PER_COVERAGE):
        '''
        ORS_List: List of computed reachable sets
        logs: Input logs to monitor
        tp: Representation type of the log
        th1: The state variable to be visualized
        T: Max time step
        vizCoverage: % of logs to be visualized
        '''
        prsr=InpParse(logFname,tp)
        logs=prsr.getLog()
        Visualize.vizMonitor(ORS_List,logs,T,th1,self.unsafeList,fname,vizCoverage)

    def vizCompMonitorLogFile(self,ORS_List_Offline,logFnameOffline,ORS_List_Online,logOnline,tp,T,th1,fname="viz_test",vizCoverage=VIZ_PER_COVERAGE):
        '''
        ORS_List_Offline: Reachable sets for the offline monitor
        logsFnameOffline: File name containing the input logs to the offline monitor
        ORS_List_Online: Reachable sets for the online monitor
        logOnline: Sythesized log by the online monitor
        tp: Representation type of the log
        th1: The state variable to be visualized
        T: Max time step
        vizCoverage: % of logs to be visualized
        '''
        prsr=InpParse(logFnameOffline,tp)
        logsOffline=prsr.getLog()
        Visualize.vizMonitorCompare(ORS_List_Offline,logsOffline,ORS_List_Online,logOnline,T,th1,self.unsafeList,fname,vizCoverage)
