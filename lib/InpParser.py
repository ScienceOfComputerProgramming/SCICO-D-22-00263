import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)
from Parameters import *

import json

from lib.GenLog import *

'''
Provides API to parse logs as input
'''

class InpParse:
    def __init__(self,logFname,tpRep='interval',tpTS='precise'):
        self.logFname=logFname # File name of the log
        self.tpRep=tpRep # representation of the log
        self.tpTS=tpTS # Type of the log
        '''
        The logs can either be `tp=interval` or `tp=zonotope`
        '''

    def getLog(self):
        if self.tpTS.lower()=='precise':
            return self.getLogPT()
        elif self.tpTS.lower()=='uncertain':
            return self.getLogUT()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Timestamp format not supported!{bcolors.ENDC}")
            exit()

    def getLogUT(self):
        if self.tpRep.lower()=='interval':
            log=self.getLogIntervalUT()
        elif self.tpRep.lower()=='zonotope':
            log=self.getLogZonoUT()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Input format not supported!{bcolors.ENDC}")
            exit()
        return log

    def getLogPT(self):
        if self.tpRep.lower()=='interval':
            log=self.getLogIntervalPT()
        elif self.tpRep.lower()=='zonotope':
            log=self.getLogZonoPT()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Input format not supported!{bcolors.ENDC}")
            exit()
        return log

    def getBehavior(self):
        if self.tpRep.lower()=='interval':
            log=self.getBehInterval()
        elif self.tpRep.lower()=='zonotope':
            log=self.getBehZono()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Input format not supported!{bcolors.ENDC}")
            exit()
        return log

    def getLogIntervalPT(self):
        '''
        Extract logs (from file) to a python datastructure.
        '''
        fl = open(DATA_PATH+self.logFname+'.mlog', "r")
        logsInterval=[]
        oldT=-1
        for ln in fl:
            logL=ln.split(':')
            if int(logL[0])<=oldT:
                print(f"{bcolors.OKCYAN}{bcolors.FAIL}Time steps have to be monotonic in the log!{bcolors.ENDC}")
                exit()
            intvl=json.loads(logL[1])
            logsInterval.append([int(logL[0]),intvl])
            oldT=int(logL[0])

        # Convert logsInterval to proper format for MoULDyS
        logs=GenLog.intvl2MlogPT(logsInterval)

        return Logs(logs,'precise','interval')

    def getLogIntervalUT(self):
        '''
        Extract logs (from file) to a python datastructure.
        '''
        fl = open(DATA_PATH+self.logFname+'.mlog', "r")
        logsInterval=[]
        oldT=-1
        for ln in fl:
            logL=ln.split(':')
            tstamp=json.loads(logL[0])

            intvl=json.loads(logL[1])
            logsInterval.append([tstamp,intvl])

        # Convert logsInterval to proper format for MoULDyS
        logs=GenLog.intvl2MlogUT(logsInterval)

        return Logs(logs,'uncertain','interval')

    def getBehInterval(self):
        '''
        Extract logs (from file) to a python datastructure.
        '''
        fl = open(DATA_PATH+self.logFname+'.mbeh', "r")
        logsInterval=[]
        for ln in fl:
            intvl=json.loads(ln)
            logsInterval.append(intvl)

        # Convert logsInterval to proper format for MoULDyS
        logs=GenLog.intvl2MBeh(logsInterval)

        return logs



    def getLogZonoUT(self):
        fl = open(DATA_PATH+self.logFname+'.mlog', "r")
        logsZono=[]
        for ln in fl:
            logL=ln.split(':')
            tstamp=json.loads(logL[0])
            if int(logL[0])<=oldT:
                print(f"{bcolors.OKCYAN}{bcolors.FAIL}Time steps have to be monotonic in the log!{bcolors.ENDC}")
                exit()
            logL2=logL[1].split(';')
            c=json.loads(logL2[0])
            G=np.array(json.loads(logL2[1]))
            logsZono.append([[tstamp[0],tstamp[1]],(c,G,[(-1,1)]*G.shape[1])])

        return Logs(logsZono,'uncertain','zonotope')

    def getLogZonoPT(self):
        fl = open(DATA_PATH+self.logFname+'.mlog', "r")
        logsZono=[]
        for ln in fl:
            logL=ln.split(':')
            logL2=logL[1].split(';')
            c=json.loads(logL2[0])
            G=np.array(json.loads(logL2[1]))
            logsZono.append([int(logL[0]),(c,G,[(-1,1)]*G.shape[1])])

        return Logs(logsZono,'precise','zonotope')

    def getBehZono(self):
        fl = open(DATA_PATH+self.logFname+'.mbeh', "r")
        logsZono=[]
        oldT=-1
        for ln in fl:
            logL2=ln.split(';')
            c=json.loads(logL2[0])
            G=np.array(json.loads(logL2[1]))
            logsZono.append([int(logL[0]),(c,G,[(-1,1)]*G.shape[1])])

        return logsZono



if False:
    InpParse('testLog.mlog','interval').getLog()
