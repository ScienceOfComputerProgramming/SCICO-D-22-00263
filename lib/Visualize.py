import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

import pickle
import matplotlib.pyplot as plt
import statistics as stat
import numpy as np
#import seaborn as sns
import pandas as pd
from matplotlib.patches import Ellipse
import sys,os
import matplotlib
from gurobipy import *
from lib.ULS_Engine.StarOperations import *
import random
import math

from Parameters import *


class Visualize:
    '''
    Visualization APIs
    '''
    def getPlotsLineFine(star,theta1,theta2):
        '''
        Returns the list of points (x,y)
        for the set star, along dimensions theta1 and theta2
        '''

        C=star[0]
        V=star[1]
        P=star[2]
        X_list=[]
        Y_list=[]

        sv=V.shape[0]
        aS=V.shape[1]

        semiDefFlag=False
        model = Model("qp")
        model.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(aS):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(aS):
            objX=objX+(predVars[i]*V[theta1][i])
        objX=C[theta1]+objX

        objY=0
        for i in range(aS):
            objY=objY+(predVars[i]*V[theta2][i])
        objY=C[theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(aS):
            a=P[i][0]
            b=P[i][1]
            if a==b:
                model.addConstr(predVars[i]==a,name)
            else:
                model.addConstr(predVars[i]>=min(a,b),name+".1")
                model.addConstr(predVars[i]<=max(a,b),name+".2")
        #-----------------------------------


        # Quadrant Specific Constraints

        # 1st Quadrant

        obj=X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(900):
            an=a/10
            if an==90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*(X),"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        #'''
        # 2nd Quadrant

        obj=X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(0,-900,-1):
            an=a/10
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 3rd Quadrant

        obj=-X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(-900,-1800,-1):
            an=a/10
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 4th Quadrant

        obj=-X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(900,1800):
            an=a/10
            if an==90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------



        #------------------------------

        #print(X_list,Y_list)
        #exit(0)
        #'''

        return (X_list,Y_list)


    def vizRS(RS_List,ORS_List,th1=0,th2=1,fname="viz_x1_x2"):
        '''
        RS_List: List of reachable sets
        ORS_List: List of overapproximate reachable sets
        '''

        RS_XY_List=[]
        ORS_XY_List=[]

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("x1")
        plt.ylabel("x2")
        #plt.xlim((-1,11))
        #plt.ylim((-1,11))

        #plt.plot([-4, 12], [0, 0], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [10, 10], color='r', linestyle='--', linewidth=1)
        #plt.plot([0, 0], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([10, 10], [-2, 12], color='r', linestyle='--', linewidth=1)
        ct=1
        period=math.floor(100/VIZ_PER_COVERAGE)

        '''for rs in ORS_List:
            if ct%period==0:
                (X,Y)=VisualizeAircraft.getPlotsLineFine(rs,th1,th2)
                ORS_XY_List.append((X,Y))
                plt.plot(X,Y,'bo',label="Perturbed",alpha=0.05)
            ct=ct+1'''

        plt.plot([-49.5, -49.5], [-10, 50],  color='r', linestyle='--', linewidth=1)
        plt.plot([11, 11], [-10, 50],  color='r', linestyle='--', linewidth=1)


        for rs in RS_List:
            if period%ct==0 or True:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                RS_XY_List.append((X,Y))
                #print(".")
                plt.plot(X,Y,'go',label="Unperturbed",alpha=0.03)
            ct=ct+1


        #plt.legend()
        plt.show()
        #plt.savefig(OUTPUT_PATH+"/"+fname+".pdf", dpi=100, bbox_inches='tight')
        plt.close()


    def vizMonitor(ORS_List,logObj,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List: List of reachable sets
        logObj: The log object
        T: Max time step
        th1: State variable to be Visualized
        unsafeList: Unsafe specification
        vizCovergae: % reachable sets to be visualized
        '''
        log=logObj.lg # The input log
        if logObj.typeTS.lower()=='precise':
            Visualize.vizMonitorPT(ORS_List,log,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE)
        elif logObj.typeTS.lower()=='uncertain':
            Visualize.vizMonitorUT(ORS_List[0],ORS_List[1],log,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE)
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Timestamp format not supported!{bcolors.ENDC}")
            exit()

    def vizMonitorPT(ORS_List,log,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List: List of reachable sets
        log: The input log object
        T: Max time step
        th1: State variable to be Visualized
        unsafeList: Unsafe specification
        vizCovergae: % reachable sets to be visualized
        '''
        time_taken=time.time()
        th2=(th1+1)%len(ORS_List[0][0])

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("State "+str(th1))


        ct=0
        period=math.floor(100/vizCovergae)
        for rs in ORS_List:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in log:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='k')

        for unS in unsafeList:
            un=unS[2]
            lb=un[th1][0]
            ub=un[th1][1]
            if lb!=-1000 and lb!=1000:
                plt.plot([-1, T], [lb, lb], color='r', linestyle='--', linewidth=1)
            elif ub!=-1000 and ub!=1000:
                plt.plot([-1, T], [ub, ub], color='r', linestyle='--', linewidth=1)

        #plt.savefig(OUTPUT_PATH+"/"+fname+"_"+id, dpi=100, bbox_inches='tight')
        plt.show()
        plt.close()

    def vizMonitorCompare(ORS_List_Offline,logOfflineObj,ORS_List_Online,logsOnlineObj,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List_Offline: Reachable sets for the offline monitor
        logsOfflineObj: Input log object to the offline monitor
        ORS_List_Online: Reachable sets for the online monitor
        logOnlineObj: Sythesized log by the online monitor
        th1: The state variable to be visualized
        T: Max time step
        vizCoverage: % of logs to be visualized
        '''
        logOffline=logOfflineObj.lg
        if logOfflineObj.typeTS=='precise':
            Visualize.vizMonitorComparePT(ORS_List_Offline,logOffline,ORS_List_Online,logsOnlineObj.lg,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE)
        elif logOfflineObj.typeTS.lower()=='uncertain':
            Visualize.vizMonitorCompareUT(ORS_List_Offline[0],ORS_List_Offline[1],logOffline,ORS_List_Online,logsOnlineObj.lg,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE)
        else:
            print(f"{bcolors.OKCYAN}{bcolors.FAIL}Timestamp format not supported!{bcolors.ENDC}")
            exit()


    def vizMonitorComparePT(ORS_List_Offline,logOffline,ORS_List_Online,logsOnline,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List_Offline: Reachable sets for the offline monitor
        logsOffline: Input logs to the offline monitor
        ORS_List_Online: Reachable sets for the online monitor
        logOnline: Sythesized log by the online monitor
        th1: The state variable to be visualized
        T: Max time step
        vizCoverage: % of logs to be visualized
        '''
        time_taken=time.time()
        th2=(th1+1)%len(ORS_List_Offline[0][0])

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("State "+str(th1))


        ct=0
        period=math.floor(100/vizCovergae)
        for rs in ORS_List_Online:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1
        for lg in logsOnline:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='k')
        ct=0
        for rs in ORS_List_Offline:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='g', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in logOffline:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='m')

        for unS in unsafeList:
            un=unS[2]
            lb=un[th1][0]
            ub=un[th1][1]
            if lb!=-1000 and lb!=1000:
                plt.plot([-1, T], [lb, lb], color='r', linestyle='--', linewidth=1)
            elif ub!=-np.inf and ub!=np.inf:
                plt.plot([-1, T], [ub, ub], color='r', linestyle='--', linewidth=1)

        #plt.savefig(OUTPUT_PATH+"/"+fname+"_"+id, dpi=100, bbox_inches='tight')
        plt.show()
        plt.close()

    def vizMonitorCompareUT(ORS_List_Offline,ORS_List_OfflineTimes,logOffline,ORS_List_Online,logsOnline,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List_Offline: Reachable sets for the offline monitor
        ORS_List_OfflineTimes: The associated uncertainties in time
        logsOffline: Input logs to the offline monitor
        ORS_List_Online: Reachable sets for the online monitor
        logsOnline: Sythesized log by the online monitor
        th1: The state variable to be visualized
        T: Max time step
        vizCoverage: % of logs to be visualized
        '''
        time_taken=time.time()
        th2=(th1+1)%len(ORS_List_Offline[0][0])

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("State "+str(th1))


        ct=0
        period=math.floor(100/vizCovergae)
        for rs in ORS_List_Online:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1
        for lg in logsOnline:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='k')


        ct=0
        for rs,t in zip(ORS_List_Offline,ORS_List_OfflineTimes):
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([t,t],[min(X),max(X)],color='g', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in logOffline:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            for lgTime in lg[0]:
                plt.plot([lgTime,lgTime],[min(X),max(X)], linestyle='-', linewidth=4,color='m')

        '''
        for rs in ORS_List_Offline:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='g', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in logOffline:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='m')
        '''


        for unS in unsafeList:
            un=unS[2]
            lb=un[th1][0]
            ub=un[th1][1]
            if lb!=-1000 and lb!=1000:
                plt.plot([-1, T], [lb, lb], color='r', linestyle='--', linewidth=1)
            elif ub!=-np.inf and ub!=np.inf:
                plt.plot([-1, T], [ub, ub], color='r', linestyle='--', linewidth=1)

        plt.savefig(OUTPUT_PATH+"/"+fname+"_"+id, dpi=100, bbox_inches='tight')
        #plt.show()
        plt.close()

    def vizMonitorUT(ORS_List,reachSetsTimes,log,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
        '''
        ORS_List: Reachable sets
        reachSetsTimes: The time uncertainties
        log: Input log
        T: Max time step
        th1: State variable to be visualized
        unsafeList: The unsafe specification
        vizCoverage: % of logs to be visualized
        '''
        time_taken=time.time()
        th2=(th1+1)%len(ORS_List[0][0])

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("State "+str(th1))


        ct=0
        period=math.floor(100/vizCovergae)
        for rs,t in zip(ORS_List,reachSetsTimes):
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([t,t],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in log:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            for lgTime in lg[0]:
                plt.plot([lgTime,lgTime],[min(X),max(X)], linestyle='-', linewidth=4,color='k')

        for unS in unsafeList:
            un=unS[2]
            lb=un[th1][0]
            ub=un[th1][1]
            if lb!=-1000 and lb!=1000:
                plt.plot([-1, T], [lb, lb], color='r', linestyle='--', linewidth=1)
            elif ub!=-1000 and ub!=1000:
                plt.plot([-1, T], [ub, ub], color='r', linestyle='--', linewidth=1)

        #plt.savefig(OUTPUT_PATH+"/"+fname+"_"+id, dpi=100, bbox_inches='tight')
        #plt.show()
        plt.savefig(OUTPUT_PATH+"/"+fname+".pdf", dpi=100, bbox_inches='tight')
        plt.close()
