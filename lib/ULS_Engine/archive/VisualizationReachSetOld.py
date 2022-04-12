'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.

Documentation: Not yet available. (TODO)
'''

import matplotlib.pyplot as plt
import numpy as np
import math
import mpmath as mp
from gurobipy import *
from operator import add

class Visualization:
    '''
    Given a generaized star, and two states; this class provides APIs to
    visualize the star.
    '''

    def __init__(self,th1,th2,st1):
        self.theta1=th1 # X Axis
        self.theta2=th2 # Y Axis
        self.star1=st1 # Star
        self.r=st1[1].shape[0]
        self.c=st1[1].shape[1]

    def getPlotsOld(self):

        C=self.star1[0]
        V=self.star1[1]
        P=self.star1[2]
        X_list=[]
        Y_list=[]

        semiDefFlag=False
        model = Model("qp")
        model.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(self.c):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(self.c):
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(self.c):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(self.c):
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

        for an in range(90):
            if an==0:
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        # 2nd Quadrant

        obj=X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(0,-90,-1):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        for an in range(-90,-180,-1):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        for an in range(90,180):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
        return (X_list,Y_list)

    def getPlotsLine(self):

        C=self.star1[0]
        V=self.star1[1]
        P=self.star1[2]
        X_list=[]
        Y_list=[]

        sv=V.shape[0]
        aS=V.shape[1]

        U=np.zeros(sv,dtype=object)

        for i in range(sv):
            s=0
            for j in range(aS):
                s=s+(V[i][j]*mp.mpi(P[j][0],P[j][1]))
            s=s+C[i]
            s_min=float(mp.nstr(s).split(',')[0][1:])
            s_max=float(mp.nstr(s).split(',')[1][:-1])
            U[i]=(s_min,s_max)

        #print(U);exit(0)

        C=self.star1[0]
        V=np.identity(sv)
        P=U

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
        for i in range(sv):
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(sv):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints

        for i in range(sv):
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

        for an in range(90):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        # 2nd Quadrant

        obj=X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for an in range(0,-90,-1):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        for an in range(-90,-180,-1):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        for an in range(90,180):
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        return (X_list,Y_list)

    def getPlotsLineFine(self):

        C=self.star1[0]
        V=self.star1[1]
        P=self.star1[2]
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
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(aS):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        return (X_list,Y_list)

    def getPlotsLineFinePred(self):

        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2][0]
        D=rs[2][1]

        sv=V.shape[0]
        aS=V.shape[1]
        n_constr=P.shape[0] # Number of Constraints on the predicate Variables

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
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(aS):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(n_constr):
            name="Pred-C"+str(i)
            cn=0
            for j in range(aS):
                name=name+"-"+str(i)
                cn=cn+(P[i][j]*predVars[j])
            a=D[i][0][0]
            b=D[i][0][1]
            model.addConstr(cn<=max(a,b),name+".1")
            model.addConstr(cn>=min(a,b),name+".2")
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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

        return (X_list,Y_list)

    def getPlotsHalfSpace(self):

        C=self.star1[0]
        V=self.star1[1]
        P=self.star1[2]
        X_list=[]
        Y_list=[]

        sv=V.shape[0]
        aS=V.shape[1]

        semiDefFlag=False
        model = Model("qp")
        model.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(sv):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(sv):
            objX=objX+(predVars[i]*V[self.theta1][i])
        objX=C[self.theta1]+objX

        objY=0
        for i in range(sv):
            objY=objY+(predVars[i]*V[self.theta2][i])
        objY=C[self.theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(sv):
            a=P[i][0]
            b=P[i][1]
            if a==b:
                model.addConstr(predVars[i]==a,name)
            else:
                model.addConstr(predVars[i]>=min(a,b),name+".1")
                model.addConstr(predVars[i]<=max(a,b),name+".2")
        #-----------------------------------

        # Rotating Half space
        '''obj=X+Y
        model.setObjective(obj,GRB.MAXIMIZE)
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
                    print('**The model cannot be solved because it is infeasible or unbounded**')
                else:
                    xVal=model.getVarByName("X").x
                    yVal=model.getVarByName("Y").x
                    print(xVal,yVal)
                    X_list.append(xVal)
                    Y_list.append(yVal)
        except:
            semiDefFlag=True
        exit(0)'''

        model.update()
        for an in range(360):
            #an=a/100
            obj=0
            if an==90 or an==270:
                obj=X
            else:
                m=math.tan(math.radians(an))
                obj=Y-(m*X)

            print(obj)
            model.setObjective(obj,GRB.MINIMIZE)
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
                        print('**The model cannot be solved because it is infeasible or unbounded**')
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
        #------------------------------------


        print(X_list,Y_list)
        return (X_list,Y_list)

    @staticmethod
    def displayPlotOld(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")

        for l in lPlots:
            X1=l[0]
            Y1=l[1]
            X2=l[2]
            Y2=l[3]
            X3=l[4]
            Y3=l[5]
            #print(X2,Y2)
            plt.plot(X1,Y1,'bo')
            plt.plot(X2,Y2,'r+')
            plt.plot(X3,Y3,'kx')

        #plt.axis('scaled')
        #plt.legend()
        #plt.savefig("Plots/"+name)
        plt.show()
        #plt.close()

    @staticmethod
    def displayPlot(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")

        for l in lPlots:
            XY=l[0]
            X1=l[1]
            Y1=l[2]
            X2=l[3]
            Y2=l[4]
            X3=l[5]
            Y3=l[6]
            #print(X2,Y2)
            for (X,Y) in XY:
                plt.plot(X,Y,'mo')
            plt.plot(X1,Y1,'bo')
            plt.plot(X2,Y2,'ro')
            plt.plot(X3,Y3,'ko')

        #plt.axis('scaled')
        #plt.legend()
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()

    @staticmethod
    def displayPlotSingle(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        X1=lPlots[0]
        Y1=lPlots[1]
        X2=lPlots[2]
        Y2=lPlots[3]


        plt.plot(X1,Y1,'r.')
        plt.plot(X2,Y2,'b.')

        #plt.axis('scaled')
        #plt.legend()
        plt.savefig("PlotsU/"+name)
        #plt.show()
        plt.close()

    def displayPlotTmp(th1,th2,lPlots):
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")
        plt.plot(lPlots[0],lPlots[1],'bo')
        #plt.axis('scaled')
        plt.legend()
        plt.show()

if False:
    C=[0,0]
    V=np.array([
    [0,-1],
    [1,0]
    ])
    P=[(-1,1), (-2,2)]
    rs=(C,V,P)
    v=Visualization(0,1,rs)
    Visualization.displayPlotTmp(0,1,v.getPlots())

if False:
    C=[0, 0, 0]
    V=np.array([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
    P=[(1,1),(1,1),(1,1)]
    rs=(C,V,P)
    C2=[0, 0, 0]
    V2=np.array([
    [1,0,0.05],
    [0,0.95,0],
    [0,0,1]
    ])
    P2=[(-0.0049999999999999975, 0.0050000000000000044), (0.0, 0.0), (0.0, 0.0)]
    rs2=(C2,V2,P2)
    v=Visualization(0,1,rs,rs2)
    v.displayPlot()

if True:
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=np.array([
    [1,1],
    [1,-1],
    ])
    D=np.array([
    [(-5,5)],
    [(-1,1)]
    ])
    rs=(C,V,(P,D))

    q=Visualization(0,1,rs)
    Visualization.displayPlotTmp(0,1,q.getPlotsLineFinePred())
