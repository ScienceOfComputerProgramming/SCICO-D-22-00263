import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
FLOW_STAR_ROOT = os.environ['FLOW_STAR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)
import numpy as np

from Parameters import *
from lib.ULS_Engine.StarOperations import *

class FlowStar:
    def __init__(self,bname):
        self.bname=bname # Benchmark Name
        if self.bname=="Anaesthesia":
            self.fname="anaesthesia3"
            self.lineNumT=7
            self.lineT='\ttime %\n'
            #self.lineNumsInitSet=[28,32]
            self.lineNumsInitSet=[28,32]
            self.lines=[1]*5
            self.lines[0]='\tcp in [% , #]\n'
            self.lines[1]='\tc1 in [% , #]\n'
            self.lines[2]='\tc2 in [% , #]\n'
            self.lines[3]='\tce in [% , #]\n'
            self.lines[4]='\tu in [% , #]\n'
        else:
            print("Not Implemented!")

    def isSafe(self,initSet,T):
        '''
        Returns reachable sets from Flow*
        '''
        #self.writeInpFile(initSet,T)
        self.runFlowStar()
        #fg=self.readOpFile()
        #return fg

    def writeInpFile(self,initSet,T):
        # Take a box hull
        boxInit=StarOp.boxHull(initSet)
        # %%%%%%%%%%%%%%%

        # Write the box-hull to the input file

        #print(os.getcwd(),PROJECT_ROOT)

        flowIpFile = open(FLOW_STAR_ROOT+'/flowstar-2.1.0/'+self.fname+'.model', 'r').readlines()

        flowIpFile[self.lineNumT]=self.lineT.replace("%",str(T*0.01))
        ct=0
        for ln in range(self.lineNumsInitSet[0],self.lineNumsInitSet[1]+1):
            flowIpFile[ln]=self.lines[ct].replace("%",str(boxInit[2][ct][0])).replace("#",str(boxInit[2][ct][1]))
            ct=ct+1

        out = open(FLOW_STAR_ROOT+'/flowstar-2.1.0/'+self.fname+'.model', 'w')
        out.writelines(flowIpFile)
        out.close()
        # %%%%%%%%%%%%%%%

    def runFlowStar(self):
        '''
        Run the flowstar code
        '''
        #print(PROJECT_ROOT+'lib/flowstar-2.1.0/flowstar-2.1.0/flowstar < '+PROJECT_ROOT+'/lib/flowstar-2.1.0/flowstar-2.1.0/'+self.fname+'.model')
        #os.system('touch '+PROJECT_ROOT+'/lib/flowstar-2.1.0/flowstar-2.1.0/outputs/'+"tester"+'.model')
        os.system(FLOW_STAR_ROOT+'/flowstar-2.1.0/flowstar < '+FLOW_STAR_ROOT+'/flowstar-2.1.0/'+self.fname+'.model')
        #os.system(FLOW_STAR_ROOT+'/flowstar-2.1.0/flowstar < '+PROJECT_ROOT+'/lib/outputs/'+self.fname+'.flow')

    def readOpFile(self):
        opFile = open(PROJECT_ROOT+'lib/counterexamples/'+self.fname+'.counterexample', 'r')
        opFileAllLines=opFile.readlines()
        opFile.close()
        #print(len(opFileAllLines))
        if len(opFileAllLines)>4:
            return 0
        else:
            return 1





if True:
    C=[0,0,0,0,0]
    V=np.array([
    [1,0,0,0,0],
    [0,1,0,0,0],
    [0,0,1,0,0],
    [0,0,0,1,0],
    [0,0,0,0,1],
    ])
    P=[(3,4),(3,4),(4,5),(3,4),(2,5)]
    #P=[(3,40),(3,40),(4,50),(3,40),(2,50)]

    initSet=(C,V,P)
    fs=FlowStar("Anaesthesia")
    fg=fs.isSafe(initSet,50)
    print("Safety: ",(fg==1))
