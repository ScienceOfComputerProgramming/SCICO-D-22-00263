import os,sys
import pickle
import time
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT)

import pickle
import matplotlib.pyplot as plt
import statistics as stat
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.patches import Ellipse
import sys,os
import matplotlib
from gurobipy import *
from lib.ULS_Engine.StarOperations import *
import random
import math

from Parameters import *


class VisualizeAircraft:
    '''
    Visualization APIs related to Aircraft Model
    '''
    def getPlotsLineFine(star,theta1,theta2):
        '''
        Returns the list of points (x,y)
        for the reachable set st1
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

    def vizX1X2(RS_List,ORS_List,fname="viz_x1_x2"):
        th1=0
        th2=1

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
        ct=0
        period=math.floor(100/VIZ_PER_COVERAGE)

        '''for rs in ORS_List:
            if ct%period==0:
                (X,Y)=VisualizeAircraft.getPlotsLineFine(rs,th1,th2)
                ORS_XY_List.append((X,Y))
                plt.plot(X,Y,'bo',label="Perturbed",alpha=0.05)
            ct=ct+1
        '''

        for rs in RS_List:
            (X,Y)=VisualizeAircraft.getPlotsLineFine(rs,th1,th2)
            RS_XY_List.append((X,Y))
            #print(".")
            plt.plot(X,Y,'go',label="Unperturbed",alpha=0.03)


        #plt.legend()
        plt.show()
        #plt.savefig(OUTPUT_PATH+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizOfflineMonitorX1(ORS_List,logs,T,fname="viz_offline_monitor_x1"):
        th1=0
        th2=1

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("x1")


        ct=0
        period=math.floor(100/VIZ_PER_COVERAGE)
        for rs in ORS_List:
            if ct%period==0:
                (X,Y)=VisualizeAircraft.getPlotsLineFine(rs,th1,th2)
                plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=4)
            ct=ct+1


        for (t,rs) in logs:
            (X,Y)=VisualizeAircraft.getPlotsLineFine(rs,th1,th2)
            plt.plot([t,t],[min(X),max(X)],color='k', linestyle='-', linewidth=4)



        plt.plot([-2, T+2],[11, 11], color='r', linestyle='--', linewidth=1)
        plt.plot([-2, T+2],[-49, -49], color='r', linestyle='--', linewidth=1)

        #plt.legend()
        plt.show()
        #plt.savefig(OUTPUT_PATH+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()
