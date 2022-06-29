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
    def __init__(self,logFname,tp='interval'):
        self.logFname=logFname # File name of the log
        self.tp=tp
        '''
        The logs can either be `tp=interval` or `tp=zonotope`
        '''

    def getLog(self):
        if self.tp.lower()=='interval':
            log=self.getLogInterval()
        elif self.tp.lower()=='zonotope':
            log=self.getLogZono()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Input format not supported!{bcolors.ENDC}")
            exit()
        return log

    def getBehavior(self):
        if self.tp.lower()=='interval':
            log=self.getBehInterval()
        elif self.tp.lower()=='zonotope':
            log=self.getBehZono()
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Input format not supported!{bcolors.ENDC}")
            exit()
        return log

    def getLogInterval(self):
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
        logs=GenLog.intvl2Mlog(logsInterval)

        return logs

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



    def getLogZono(self):
        fl = open(DATA_PATH+self.logFname+'.mlog', "r")
        logsZono=[]
        oldT=-1
        for ln in fl:
            logL=ln.split(':')
            if int(logL[0])<=oldT:
                print(f"{bcolors.OKCYAN}{bcolors.FAIL}Time steps have to be monotonic in the log!{bcolors.ENDC}")
                exit()
            logL2=logL[1].split(';')
            c=json.loads(logL2[0])
            G=np.array(json.loads(logL2[1]))
            logsZono.append([int(logL[0]),(c,G,[(-1,1)]*G.shape[1])])

        return logsZono

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
