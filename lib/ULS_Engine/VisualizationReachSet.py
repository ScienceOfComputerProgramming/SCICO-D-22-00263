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
import io
from PIL import Image

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

    def getPlotsLineFineList(RS_list,s1,s2):
        pltList=[]
        for rs in RS_list:
            (X,Y)=Visualization(s1,s2,rs).getPlotsLineFine()
            pltList.append((X,Y))
        return pltList

    def getPlotsLineFinePred(self):

        C=self.star1[0] # The center is always assumed to be 0 as of now
        V=self.star1[1]
        P=self.star1[2][0]
        D=self.star1[2][1]

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
                        0+0
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
                        0+0
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
                        0+0
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
                        0+0
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
        an=-1
        m=math.tan(math.radians(an))
        obj=Y-(m*X)
        #obj=X+Y
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
        Visualization.displayPlotTmp(0,1,(X_list,Y_list))
        exit(0)

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
        #plt.savefig("Plots/"+name)
        plt.show()
        plt.close()

    @staticmethod
    def displayPlotORIGINAL(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")

        for l in lPlots:
            XY=l[0] # Random Samples
            X2=l[1] # RS X
            Y2=l[2] # RS Y
            X3=l[3] # ORS X
            Y3=l[4] # ORS Y
            XY2=l[5] # Random Predicates
            #print(X2,Y2)
            for (X,Y) in XY:
                plt.plot(X,Y,'mo')
            for (X,Y) in XY2:
                plt.plot(X,Y,'co',alpha=0.02)
            plt.plot(X2,Y2,'bo')
            plt.plot(X3,Y3,'k.')

        #plt.axis('scaled')
        #plt.legend()
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()

    @staticmethod
    def displayPlot(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")

        for l in lPlots:
            X1=l[0] # RS X
            Y1=l[1] # RS Y
            X2=l[2] # ORS_compact X
            Y2=l[3] # ORS_compact Y
            X3=l[4] # ORS X
            Y3=l[5] # ORS Y
            #print(X2,Y2)
            plt.plot(X1,Y1,'bo')
            plt.plot(X3,Y3,'ko')
            plt.plot(X2,Y2,'ro',alpha=0.04)

        #plt.axis('scaled')
        #plt.legend()
        #plt.savefig("Plots/"+name)
        plt.show()
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
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()

    @staticmethod
    def displayPlotMultIS(th1,th2,lPlots,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        X1=lPlots[0]
        Y1=lPlots[1]
        X2=lPlots[2]
        Y2=lPlots[3]
        X3=lPlots[4]
        Y3=lPlots[5]
        X4=lPlots[6]
        Y4=lPlots[7]


        plt.plot(X1,Y1,'ro')
        plt.plot(X2,Y2,'bo')
        plt.plot(X3,Y3,'ko')
        plt.plot(X4,Y4,'co')

        #plt.axis('scaled')
        #plt.legend()
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()

    def displayPlotTmp(th1,th2,lPlots,name):
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        #plt.plot(X2,Y2,'r+',label="Perturbed")
        plt.plot(lPlots[0],lPlots[1],'bo')
        #plt.axis('scaled')
        plt.legend()
        plt.show()

    def displayPlotList(th1,th2,lPlots,XY,name):
        plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        plt.plot(lPlots[0][0],lPlots[0][1],'bo')
        plt.plot(lPlots[1][0],lPlots[1][1],'co')
        plt.plot(lPlots[2][0],lPlots[2][1],'mo')
        plt.plot(lPlots[3][0],lPlots[3][1],'ko')
        plt.plot(XY[0],XY[1],'ro',alpha=0.1)
        #plt.plot(lPlots[3][0],lPlots[3][1],'bo')
        #plt.legend()
        #plt.show()
        plt.savefig("Plots/"+name)
        plt.close()

    def getPlotList(th1,th2,lPlots,XY,name):
        #plt.axes()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        plt.plot(lPlots[0][0],lPlots[0][1],'bo')
        plt.plot(lPlots[1][0],lPlots[1][1],'co')
        plt.plot(lPlots[2][0],lPlots[2][1],'mo')
        plt.plot(lPlots[3][0],lPlots[3][1],'ko')
        plt.plot(XY[0],XY[1],'ro',alpha=0.1)
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        im = Image.open(buf)
        plt.savefig("Plots/"+name)
        plt.close()
        #im.show()
        return im

    def getPlotAll(th1,th2,lPlots,flow,name):
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        XY=lPlots[0]
        X1=lPlots[1][0]
        Y1=lPlots[1][1]
        X2=lPlots[1][2]
        Y2=lPlots[1][3]
        X3=lPlots[1][4]
        Y3=lPlots[1][5]
        X4=lPlots[1][6]
        Y4=lPlots[1][7]
        for (X,Y) in XY:
            plt.plot(X,Y,'mo')
        plt.plot(X1,Y1,'bo') # Nominal
        plt.plot(X2,Y2,'co') # Zono reduction
        plt.plot(X3,Y3,'go') # Interval reduction
        plt.plot(X4,Y4,'ro',alpha=0.4) # No reduction

        plt.plot(flow[0],flow[1],color='black',linewidth=6) # Flow*

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        im = Image.open(buf)
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()
        return im

    def getPlotTwo(th1,th2,lPlots,name):
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        XY=lPlots[0]
        X1=lPlots[1][0]
        Y1=lPlots[1][1]
        X2=lPlots[1][2]
        Y2=lPlots[1][3]
        X3=lPlots[1][4]
        Y3=lPlots[1][5]

        for (X,Y) in XY:
            plt.plot(X,Y,'mo')
        plt.plot(X1,Y1,'bo')
        plt.plot(X2,Y2,'ro')
        plt.plot(X3,Y3,'go')

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        im = Image.open(buf)
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()
        return im

    def getPlotOrdComp(th1,th2,lPlots,name):
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        X1=lPlots[0]
        Y1=lPlots[1]
        XY=lPlots[2] #magenta
        XY2=lPlots[3] #green

        for (X,Y) in XY:
            plt.plot(X,Y,'mo')

        for (X,Y) in XY2:
            plt.plot(X,Y,'go')

        plt.plot(X1,Y1,'bo')


        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        im = Image.open(buf)
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()
        return im

    def getPlotPCA(th1,th2,lPlots,name):
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))

        XY=lPlots[0]
        XY2=lPlots[1]
        X1=lPlots[2][0]
        Y1=lPlots[2][1]
        X2=lPlots[2][2]
        Y2=lPlots[2][3]
        X3=lPlots[2][4]
        Y3=lPlots[2][5]

        for (X,Y) in XY2:
            plt.plot(X,Y,'mo')

        plt.plot(X1,Y1,'bo')
        plt.plot(X2,Y2,'ro')
        plt.plot(X3,Y3,'co')

        for (X,Y) in XY:
            plt.plot(X,Y,'ko')

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        im = Image.open(buf)
        plt.savefig("Plots/"+name)
        #plt.show()
        plt.close()
        return im

    @staticmethod
    def displayPlotEgHeuristics(th1,th2,ac,intvl,zono):
        #plt.axes()
        (X1,Y1)=Visualization(th1,th2,ac).getPlotsLineFine()
        (X2,Y2)=Visualization(th1,th2,intvl).getPlotsLineFine()
        (X3,Y3)=Visualization(th1,th2,zono).getPlotsLineFine()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))


        plt.plot(X1,Y1,'ro') # Actual
        plt.plot(X2,Y2,'bo') # Interval
        plt.plot(X3,Y3,'go') # Zono

        #plt.axis('scaled')
        #plt.legend()
        #plt.savefig("Plots/"+name)
        plt.show()
        plt.close()

    def displayRSnUnsafe(RS,unsafe,th1=0,th2=1,name="RS"):
        (RS_X,RS_Y)=Visualization(th1,th2,RS).getPlotsLineFine()
        (U_X,U_Y)=Visualization(th1,th2,unsafe).getPlotsLineFine()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        plt.plot(RS_X,RS_Y,'bo')
        plt.plot(U_X,U_Y,'ro')
        #plt.savefig("Plots/"+name)
        plt.show()
        plt.close()

    def displayCompHeu(RS_lin,RS_har,RS_eq,th1=0,th2=1,name="RS"):
        (RS_lin_X,RS_lin_Y)=Visualization(th1,th2,RS_lin).getPlotsLineFine()
        (RS_har_X,RS_har_Y)=Visualization(th1,th2,RS_har).getPlotsLineFine()
        (RS_eq_X,RS_eq_Y)=Visualization(th1,th2,RS_eq).getPlotsLineFine()
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        plt.plot(RS_lin_X,RS_lin_Y,'bo')
        plt.plot(RS_har_X,RS_har_Y,'go')
        plt.plot(RS_eq_X,RS_eq_Y,'ko')
        #plt.savefig("Plots/"+name)
        plt.show()
        plt.close()

if False:
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=np.array([
    [1,1],
    [1,-1],
    [-4,4],
    ])
    D=np.array([
    [(-1,1)],
    [(-1,1)],
    [(-1,1)]
    ])
    rs=(C,V,(P,D))

    q=Visualization(0,1,rs)
    Visualization.displayPlotTmp(0,1,q.getPlotsLineFinePred())

if False:
    C=[100,100]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=[(-1,1),(-1,1)]
    Visualization(0,1,(C,V,P)).getPlotsHalfSpace()
