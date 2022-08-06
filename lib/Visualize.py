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


    def vizTest(ORS_List,log,T,fname="viz_test"):
        time_taken=time.time()
        th1=0
        th2=1

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("h")


        ct=0
        period=math.floor(100/VIZ_PER_COVERAGE)
        for rs in ORS_List:
            if ct%period==0:
                (X,Y)=Visualize.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1

        for lg in log:
            (X,Y)=Visualize.getPlotsLineFine(lg[1],th1,th2)
            plt.plot([lg[0],lg[0]],[min(X),max(X)], linestyle='-', linewidth=4,color='k')

        plt.plot([-1, T], [0.5, 0.5], color='r', linestyle='--', linewidth=1)
        idT=str(time.time()).split('.')
        id=''
        PR=4
        for i in range(PR-1,-1,-1):
            id=id+idT[0][-i]
        id=id+str('-')+idT[1][0]+idT[1][1]

        #plt.savefig(OUTPUT_PATH+"/"+fname+"_"+id, dpi=100, bbox_inches='tight')
        plt.close()

    def vizMonitor(ORS_List,log,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
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

    def vizMonitorCompare(ORS_List_Offline,logOffline,ORS_List_Online,logsOnline,T,th1,unsafeList,fname="viz_test",vizCovergae=VIZ_PER_COVERAGE):
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
