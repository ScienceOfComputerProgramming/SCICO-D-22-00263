import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

from Parameters import *
import random
from lib.ULS_Engine.ComputeU import *
from lib.ULS_Engine.StarOperations import *
import copy
import json

class Logs:
    def __init__(self,lg,typeTS="precise",typeRep="zono"):
        self.lg=lg
        self.typeTS=typeTS
        self.typeRep=typeRep

    def getTimestampType(self):
        return self.typeTS

    def getRepresentationType(self):
        return self.typeRep

    def getLog(self):
        return self.lg

class GenLog:
    '''
    Given a linear uncertain system and an initial set, generate logs for monitoring
    '''
    def __init__(self,A,Er,initialSet,T,pr=PROBABILITY_LOG,delta=0):
        self.A=A # The nominal dynamics matrix
        self.Er=Er # Uncertainties
        self.initialSet=initialSet # Initial Set
        self.T=T # The maximum time horizon
        self.pr=pr # The probability of logging at a given step is p/100
        self.delta=delta # Uncertaity in time

    def getAccLog(self):
        '''
        Generate a log:
            - Take a single point from initial set.
            - Assume a fixed valuation of the uncertainties.
            - At each step, log values with probability p.

        Output:
            1. The generated log (a two tuple of time and set). Note that
            this returns an accurate log, with no uncertainties.
            2. The actual reachable set at all time steps.
        '''
        log=[] # log: a list of reachable points.

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

        return (Logs(log,'precise','zono'),actualBehavior)

    def getLog(self):
        if self.delta==0:
            return self.getLogPTime()
        else:
            return self.getLogUTime()


    def getLogPTime(self):
        '''
        Get a log with uncertaity in it:
            - Generate an accurate log.
            - Add uncertaity to the reachable set.

        Note: The purpose of this API is to emulate the actual behavior
        for the online monitoring algorithms logger function.
        That is, when the logger function in online monitoring is invoked,
        the actual logging is done from this simulated behavior.

        Output:
            1. The generated log (a two tuple of time and set). Note that
            this returns an accurate log, with uncertainties.
            2. The actual reachable set at all time steps.
        '''
        print(">>STATUS: Generating a log  . . .")
        time_taken=time.time()
        accLogObj,actualBehavior=self.getAccLog()
        accLog=accLogObj.lg
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

        return (Logs(log,'precise','zono'),actualBehaviorUn)

    def getLogUTime(self):
        '''
        Get a log with uncertaity in it:
            - Generate an accurate log.
            - Add uncertaity to the reachable set.
            - Add uncertaity to timestamps

        Note: The purpose of this API is to emulate the actual behavior
        for the online monitoring algorithms logger function.
        That is, when the logger function in online monitoring is invoked,
        the actual logging is done from this simulated behavior.

        Output:
            1. The generated log (a two tuple of time and set). Note that
            this returns an accurate log, with uncertainties.
            2. The actual reachable set at all time steps.
        '''
        print(">>STATUS: Generating a log  . . .")
        time_taken=time.time()
        accLogObj,actualBehavior=self.getAccLog()
        accLog=accLogObj.lg
        log=[]
        actualBehaviorUn=[]

        C=[0]*self.A.shape[0]
        V=np.identity(self.A.shape[0])
        P=[(-0.4,0.4)]*self.A.shape[0]
        unitBox=(C,V,P)

        lenLog=len(accLog)
        utime=self.delta

        fg=0
        for i in range(lenLog):
            lg=accLog[i]
            if fg!=0:
                rsLg=CompU.addStars(lg[1],unitBox) # Additing uncertaity to reachable set

                # Need to add uncertaity in time
                amtUTP=random.randint(0,utime/2)
                amtUTN=random.randint(0,utime/2)
                if i>amtUTN and accLog[i-1][0]>=lg[0]-amtUTN:
                    amtUTN=lg[0]-accLog[i-1][0]-1
                if i<lenLog-1 and accLog[i+1][0]<=lg[0]+amtUTP:
                    amtUTP=accLog[i+1][0]-lg[0]-1
                log.append([[lg[0]-amtUTN,lg[0]+amtUTP],rsLg])

            else:
                log.append([[lg[0]],lg[1]])
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

        return (Logs(log,'uncertain','zono'),actualBehaviorUn)


    def getLogs(self,N=5):
        '''
        Note: The purpose of this API is to emulate the actual behavior
        for the online monitoring algorithms logger function.
        That is, when the logger function in online monitoring is invoked,
        the actual logging is done from this simulated behavior.

        Output:
            1. The generated log (a two tuple of time and set). Note that
            this returns an accurate log, with uncertainties.
            2. The actual reachable set at all time steps.
        '''
        logs=[]
        print(">>STATUS: Generating logs  . . .")
        time_taken=time.time()
        for i in range(N):
            logs.append(self.getLog())
        print("\t>>STATUS: Time Taken: ",time.time()-time_taken)
        print(">>STATUS: Logs generated!")

        return logs

    def intvl2MlogPT(intvlLog):
        '''
        Convert interval logs to proper format
        '''
        if len(intvlLog)==0:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Empty log!{bcolors.ENDC}")
            exit()
        C=[0]*len(intvlLog[0][1])
        V=np.identity(len(intvlLog[0][1]))

        logs=[]
        for itvl in intvlLog:
            P=[]
            for p in itvl[1]:
                P.append((p[0],p[1]))
            rs=(copy.copy(C),copy.copy(V),copy.copy(P))
            logs.append([copy.copy(itvl[0]),copy.copy(rs)])
        return logs

    def intvl2MlogUT(intvlLog):
        '''
        Convert interval logs to proper format
        '''
        if len(intvlLog)==0:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Empty log!{bcolors.ENDC}")
            exit()
        C=[0]*len(intvlLog[0][1])
        V=np.identity(len(intvlLog[0][1]))

        logs=[]
        for itvl in intvlLog:
            P=[]
            for p in itvl[1]:
                P.append((p[0],p[1]))
            rs=(copy.copy(C),copy.copy(V),copy.copy(P))
            logs.append([[copy.copy(itvl[0][0]),copy.copy(itvl[0][1])],copy.copy(rs)])
        return logs

    def intvl2MBeh(intvlLog):
        '''
        Convert interval logs to proper format
        '''
        if len(intvlLog)==0:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Empty log!{bcolors.ENDC}")
            exit()
        C=[0]*len(intvlLog[0])
        V=np.identity(len(intvlLog[0]))

        logs=[]
        for itvl in intvlLog:
            P=[]
            for p in itvl:
                P.append((p[0],p[1]))
            rs=(copy.copy(C),copy.copy(V),copy.copy(P))
            logs.append(copy.copy(rs))
        return logs

    def getLogFile(self,fname,tp='interval'):
        if self.delta==0:
            return self.getLogFilePT(fname,tp)
        else:
            return self.getLogFileUT(fname,tp)

    def getLogFileUT(self,fname,tp='interval'):
        (log,actualBehaviorUn)=self.getLog()
        if tp.lower()=='interval':
            GenLog.log2FileIntUT(log.lg,fname)
            GenLog.behavior2FileInt(actualBehaviorUn,fname)
        elif tp.lower()=='zonotope':
            GenLog.log2FileZonoUT(log.lg,fname)
            GenLog.behavior2FileZono(actualBehaviorUn,fname)
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Cannot generate logs of this type!{bcolors.ENDC}")
        return (log,actualBehaviorUn)

    def getLogFilePT(self,fname,tp='interval'):
        (log,actualBehaviorUn)=self.getLog()
        if tp.lower()=='interval':
            GenLog.log2FileIntPT(log.lg,fname)
            GenLog.behavior2FileInt(actualBehaviorUn,fname)
        elif tp.lower()=='zonotope':
            GenLog.log2FileZonoPT(log.lg,fname)
            GenLog.behavior2FileZono(actualBehaviorUn,fname)
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Cannot generate logs of this type!{bcolors.ENDC}")
        return (log,actualBehaviorUn)


    def log2FileIntPT(logs,fname):
        lines=[]
        for lg in logs:
            tstamp=lg[0]
            rs=lg[1]
            box=StarOp.boxHull(rs)
            P=[]
            for p in box[2]:
                P.append([p[0],p[1]])
            ln=str(tstamp)+": "+json.dumps(P)
            lines.append(ln)
        with open(DATA_PATH+'/'+fname+'.mlog', 'w') as f:
            for item in lines:
                f.write("%s\n" % item)

    def log2FileIntUT(logs,fname):
        lines=[]
        for lg in logs:
            tstamp=lg[0]
            if len(tstamp)>1:
                tstamp_lb=lg[0][0]
                tstamp_ub=lg[0][1]
            else:
                tstamp_lb=lg[0][0]
                tstamp_ub=lg[0][0]
            rs=lg[1]
            box=StarOp.boxHull(rs)
            P=[]
            for p in box[2]:
                P.append([p[0],p[1]])
            ln="["+str(tstamp_lb)+","+str(tstamp_ub)+"] : "+json.dumps(P)
            lines.append(ln)
        with open(DATA_PATH+'/'+fname+'.mlog', 'w') as f:
            for item in lines:
                f.write("%s\n" % item)

    def behavior2FileInt(actualBehaviorUn,fname):
        acBeh=[]
        for beh in actualBehaviorUn:
            box=StarOp.boxHull(beh)
            P=[]
            for p in box[2]:
                P.append([p[0],p[1]])
            ln=json.dumps(P)
            acBeh.append(ln)
        with open(DATA_PATH+'/'+fname+'.mbeh', 'w') as f:
            for item in acBeh:
                f.write("%s\n" % item)

    def log2FileZonoPT(logs,fname):
        lines=[]
        for lg in logs:
            tstamp=lg[0]
            rs=lg[1]
            Z=StarOp.star2Zono(rs)
            c=Z[0].reshape(-1).tolist()
            G=Z[1].tolist()
            ln=str(tstamp)+": "+json.dumps(c)+"; "+json.dumps(G)
            lines.append(ln)
        with open(DATA_PATH+'/'+fname+'.mlog', 'w') as f:
            for item in lines:
                f.write("%s\n" % item)

    def log2FileZonoUT(logs,fname):
        lines=[]
        for lg in logs:
            tstamp=lg[0]
            if len(tstamp)>1:
                tstamp_lb=lg[0][0]
                tstamp_ub=lg[0][1]
            else:
                tstamp_lb=lg[0][0]
                tstamp_ub=lg[0][0]
            rs=lg[1]
            Z=StarOp.star2Zono(rs)
            c=Z[0].reshape(-1).tolist()
            G=Z[1].tolist()
            ln="["+str(tstamp_lb)+","+str(tstamp_ub)+"] : "+json.dumps(c)+"; "+json.dumps(G)
            lines.append(ln)
        with open(DATA_PATH+'/'+fname+'.mlog', 'w') as f:
            for item in lines:
                f.write("%s\n" % item)

    def behavior2FileZono(actualBehaviorUn,fname):
        lines=[]
        for beh in actualBehaviorUn:
            Z=StarOp.star2Zono(beh)
            c=Z[0].reshape(-1).tolist()
            G=Z[1].tolist()
            ln=json.dumps(c)+"; "+json.dumps(G)
            lines.append(ln)
        with open(DATA_PATH+'/'+fname+'.mbeh', 'w') as f:
            for item in lines:
                f.write("%s\n" % item)
