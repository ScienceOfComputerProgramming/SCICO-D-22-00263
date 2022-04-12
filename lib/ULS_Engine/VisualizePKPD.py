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
from StarOperations import *

from Parameters import *


class VisualizePKPD:
    '''
    Visualization APIs related to PKPD Model
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

    def vizC1C2(RS_List,ORS_List,fname="viz_C1_C2"):
        th1=1
        th2=2

        RS_XY_List=[]
        ORS_XY_List=[]

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("C1")
        plt.ylabel("C2")
        plt.xlim((-1,11))
        plt.ylim((-1,11))

        plt.plot([-4, 12], [0, 0], color='r', linestyle='--', linewidth=1)
        plt.plot([-4, 12], [10, 10], color='r', linestyle='--', linewidth=1)
        plt.plot([0, 0], [-2, 12], color='r', linestyle='--', linewidth=1)
        plt.plot([10, 10], [-2, 12], color='r', linestyle='--', linewidth=1)

        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            ORS_XY_List.append((X,Y))
            plt.plot(X,Y,'bo',label="Perturbed",alpha=0.05)

        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            RS_XY_List.append((X,Y))
            #print(".")
            plt.plot(X,Y,'go',label="Unperturbed",alpha=0.03)


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizCpCe(RS_List,ORS_List,fname="viz_Cp_Ce"):
        th1=0
        th2=3

        RS_XY_List=[]
        ORS_XY_List=[]

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Cp")
        plt.ylabel("Ce")
        #plt.xlim((-1,11))
        #plt.ylim((-1,11))

        #plt.plot([-4, 12], [1, 1], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [8, 8], color='r', linestyle='--', linewidth=1)
        #plt.plot([1, 1], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([6, 6], [-2, 12], color='r', linestyle='--', linewidth=1)

        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            ORS_XY_List.append((X,Y))
            plt.plot(X,Y,'bo',label="Perturbed")

        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            RS_XY_List.append((X,Y))
            #print(".")
            plt.plot(X,Y,'go',label="Unperturbed")


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizCp_old(RS_List,ORS_List,fname="viz_Cp"):
        th1=0
        th2=3

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Cp")
        #plt.xlim((-1,11))
        #plt.ylim((-1,11))

        #plt.plot([-4, 12], [1, 1], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [8, 8], color='r', linestyle='--', linewidth=1)
        #plt.plot([1, 1], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([6, 6], [-2, 12], color='r', linestyle='--', linewidth=1)

        ct=0
        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct]*len(X),X,'bo',label="Perturbed")
            ct=ct+1

        ct=0
        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct]*len(X),X,'go',label="Unperturbed")
            ct=ct+1


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizCp(RS_List,ORS_List,fname="viz_Cp"):
        th1=0
        th2=3

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Cp")
        plt.xlim((-1,21))
        plt.ylim((0,7))

        #plt.plot([1, 1], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([6, 6], [-2, 12], color='r', linestyle='--', linewidth=1)

        ct=0
        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=10)
            ct=ct+1

        ct=0
        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(X),max(X)],color='cyan', linestyle='-', linewidth=4)
            ct=ct+1

        plt.plot([-4, 25], [1, 1], color='r', linestyle='--', linewidth=1)
        plt.plot([-4, 25], [6, 6], color='r', linestyle='--', linewidth=1)

        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizCe_old(RS_List,ORS_List,fname="viz_Ce"):
        th1=0
        th2=3

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Ce")
        #plt.xlim((-1,11))
        #plt.ylim((-1,11))

        #plt.plot([-4, 12], [1, 1], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [8, 8], color='r', linestyle='--', linewidth=1)
        #plt.plot([1, 1], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([6, 6], [-2, 12], color='r', linestyle='--', linewidth=1)

        ct=0
        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct]*len(Y),Y,'bo',label="Perturbed")
            ct=ct+1

        ct=0
        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct]*len(Y),Y,'go',label="Unperturbed")
            ct=ct+1


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizCe(RS_List,ORS_List,fname="viz_Ce"):
        th1=0
        th2=3

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("Ce")
        plt.xlim((-1,21))
        plt.ylim((0,9))

        #plt.plot([-4, 12], [1, 1], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [8, 8], color='r', linestyle='--', linewidth=1)
        #plt.plot([1, 1], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([6, 6], [-2, 12], color='r', linestyle='--', linewidth=1)

        ct=0
        for rs in ORS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(Y),max(Y)],color='b', linestyle='-', linewidth=10)
            ct=ct+1


        ct=0
        for rs in RS_List:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(Y),max(Y)],color='cyan', linestyle='-', linewidth=4)
            ct=ct+1

        plt.plot([-4, 25], [1, 1], color='r', linestyle='--', linewidth=1)
        plt.plot([-4, 25], [8, 8], color='r', linestyle='--', linewidth=1)


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()

    def vizComp(RS_Top,RS_bot,fname="viz_compare"):
        th1=1
        th2=2

        RS_XY_List=[]
        ORS_XY_List=[]

        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("Time")
        plt.ylabel("C1")
        plt.xlim((-1,21))
        plt.ylim((-2,12))

        plt.plot([-4, 25], [0, 0], color='r', linestyle='--', linewidth=1)
        plt.plot([-4, 25], [10, 10], color='r', linestyle='--', linewidth=1)
        #plt.plot([-4, 12], [10, 10], color='r', linestyle='--', linewidth=1)
        #plt.plot([0, 0], [-2, 12], color='r', linestyle='--', linewidth=1)
        #plt.plot([10, 10], [-2, 12], color='r', linestyle='--', linewidth=1)

        ct=0
        for rs in RS_Top:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(X),max(X)],color='b', linestyle='-', linewidth=10)
            ct=ct+1

        ct=0
        for rs in RS_bot:
            (X,Y)=VisualizePKPD.getPlotsLineFine(rs,th1,th2)
            plt.plot([ct,ct],[min(X),max(X)],color='cyan', linestyle='-', linewidth=4)
            ct=ct+1


        #plt.legend()
        plt.savefig(PKPD_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
        plt.close()
