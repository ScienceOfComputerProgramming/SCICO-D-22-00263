'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com
- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.
Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import mpmath as mp
from gurobipy import *
import sys
import time
import math
import random

#from pypolycontain.lib.zonotope import zonotope,zonotope_order_reduction_outer,zonotope_distance_cocenter
from PIL import Image

from lib.ULS_Engine.VisualizationReachSet import *
from lib.ULS_Engine.ComputeU import *
from lib.ULS_Engine.SamplingMet import *
from lib.ULS_Engine.Gridding import *

import lib.ULS_Engine.Profiling

BIGM=0.001
EPSILON=1e-10
PRED_EP=1e-3
INTERVAL=1000000
RED_INT_ZONO=35
RED_INT_INTRVL=35
SAMPLES=20
VLENGTH=17


class Split:
    '''
    Computes the reachable set of a given linear discrete dynamical system A
    form the initial set Theta, upto time T
    '''

    def __init__(self,A,E,theta,T):
        self.A=A # Linear system (without perturbation)
        self.Er=E
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] has a perturbation of -(1-a)% to (b-1)%
        For eg: [0.9,1,1] means a perturbation of -10% and +10%
        '''
        self.Theta=theta # The initial set
        self.T=T #Max number of steps
        self.n=A.shape[0]
        self.Ac=self.computeCenter()
        self.methodName="Split"

    def computeCenter(self):
        '''
        Computes the center point matrix of the interval matrix self.A
        '''
        Ac=np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                if (i,j) in self.Er:
                    a=float(self.A[i][j]*self.Er[(i,j)][0])
                    b=float(self.A[i][j]*self.Er[(i,j)][1])
                    c=(a+b)/2
                    Ac[i][j]=c
                else:
                    Ac[i][j]=self.A[i][j]
        return Ac

    def computeU_Interval(self,rs):
        '''
        Computes the effect of uncertainty on the reachable set
        using Interval arithmetic
        '''
        A_tilde=self.computeUncertainMat()
        diff=A_tilde-self.Ac
        U=np.matmul(diff,rs)
        return U

    def computeUncertainMat(self):
        '''
        Computes the interval uncertain matrix
        '''
        A_tilde=np.zeros((self.n,self.n),dtype=object)
        for i in range(self.n):
            for j in range(self.n):
                A_tilde[i][j]=self.A[i][j]
        for key in self.Er:
            a=float(self.Er[key][0]*self.A[key[0]][key[1]])
            b=float(self.Er[key][1]*self.A[key[0]][key[1]])
            #print(a,b)
            #print(mp.mpi(float(a),float(b)))
            #A_tilde[key[0]][key[1]]=mp.mpi(float(self.Er[key][0]*self.A[key[0]][key[1]]),float(self.Er[key][1]*self.A[key[0]][key[1]]))
            #A_tilde[key[0]][key[1]]=mp.mpi(self.Er[key][0]*self.A[key[0]][key[1]],self.Er[key][1]*self.A[key[0]][key[1]])
            A_tilde[key[0]][key[1]]=mp.mpi(min(a,b),max(a,b))

        #print(A_tilde)
        return A_tilde

    def computeU(self,rs):
        '''
        Computes the effect of uncertainty on the reachable set
        using optimization.
        '''

        semiDefFlag=False

        model = Model("qp")
        model.setParam( 'OutputFlag', False )
        #model.params.Presolve=0

        # Create Perturbation Variables
        faultVars=[]
        for key in self.Er:
            name="Pert"+str(key)
            faultVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create Reachable Set Variables
        reachVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            reachVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Add the Perturbation Constraints


        model.optimize() # Updating the model, to use some of the functions

        for var in faultVars:
            kN=var.varName
            i=int(kN.split(',')[0][5:])
            j=int(kN.split(',')[1][:-1])
            k=(i,j)
            name="Pert-C"+str(k)
            model.addConstr(var>=self.Er[k][0],name+".1")
            model.addConstr(var<=self.Er[k][1],name+".2")
        #---------------------------------

        # Add Initial Set Constraints
        for i in range(self.n):
            name="ReachSet-C-"+str(i)
            if (isinstance(rs[i][0],int)) or (isinstance(rs[i][0],np.int64)) or (isinstance(rs[i][0],float)) or (isinstance(rs[i][0],np.float128)):
                model.addConstr(reachVars[i]==rs[i][0],name)
            else:
                a=float(mp.nstr(rs[i][0]).split(',')[0][1:])
                b=float(mp.nstr(rs[i][0]).split(',')[1][:-1])
                if a==b:
                    model.addConstr(reachVars[i]==a,name)
                else:
                    model.addConstr(reachVars[i]>=min(a,b),name+".1")
                    model.addConstr(reachVars[i]<=max(a,b),name+".2")
        #---------------------------------

        # Prepare the objective functions and get min max

        U=np.zeros((self.n,1),dtype=object)

        for i in range(self.n):
            obj=0
            for j in range(self.n):
                if (i,j) in self.Er:
                    pertV=model.getVarByName("Pert"+str((i,j)))
                    obj=obj+(((self.A[i][j]*pertV)-self.Ac[i][j])*reachVars[j])
                else:
                    obj=obj+((self.A[i][j]-self.Ac[i][j]))*reachVars[j]

            # Obtain Minimum
            mn=-9890
            #print(obj," (Min)\n")
            model.setObjective(obj,GRB.MINIMIZE)
            #model.params.Presolve=0
            #model.write("logMin.lp")
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
                        mn=obj.getValue()
            except:
                semiDefFlag=True


            #print("Min: ",mn)
            #-------------------------------

            # Obtain Maximum
            #print(obj," (Max)\n")
            mx=9890
            model.setObjective(obj,GRB.MAXIMIZE)
            #model.params.Presolve=0
            #model.write("logMax.lp")
            try:
                model.optimize()
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        for v in model.getVars():
                            print('%s %g' % (v.varName, v.x))
                        print('Obj: %g' % obj.getValue())
            except:
                #print("Caught it!!")
                semiDefFlag=True


            #print("Max: ",mx)
            #----------------------------------

            #U[i][0]=mp.mpi(mn,mx)
            U[i][0]=(mn,mx)
            #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ",i)

        if semiDefFlag==True:
            UInterval=self.computeU_Interval(rs)
            for i in range(self.n):
                lb=U[i][0][0]
                ub=U[i][0][1]
                if lb==-9890:
                    lb=np.float(mp.nstr(UInterval[i][0]).split(',')[0][1:])
                if ub==9890:
                    ub=np.float(mp.nstr(UInterval[i][0]).split(',')[1][:-1])
                #print(lb,ub)
                if lb!=-9890:
                    a=np.float(mp.nstr(UInterval[i][0]).split(',')[0][1:])
                    if lb<a:
                        print("Really?! Min Error!!")
                        print("Value: ",obj.getValue())
                        print("Interval: ",a)
                        '''print("i: ",i)
                        for j in range(self.n):
                            if (i,j) in self.Er:
                                pertV=model.getVarByName("Pert"+str((i,j)))
                                obj=obj+((self.A[i][j]*pertV)-self.Ac[i][j])
                            else:
                                obj=obj+(self.A[i][j]-self.Ac[i][j])
                        obj=obj*reachVars[i]
                        print(obj)
                        model.setObjective(obj,GRB.MINIMIZE)
                        model.optimize()
                        print("Value: ",obj.getValue())
                        print("Interval: ",a,UInterval[i][0])
                        print(rs)
                        print("-")
                        print(self.A[i])
                        print("-")
                        print(self.Ac[i])
                        print("-")
                        print(self.Er)
                        print("-")
                        A_tilde=self.computeUncertainMat()
                        print(A_tilde[i])
                        print("-")
                        print(self.Ac[i])
                        diff=A_tilde-self.Ac
                        print("-")
                        print(diff[i])
                        UW=np.matmul(diff,rs)
                        print(UW[i])
                        model.write('dump.lp')'''
                        #exit(0)
                if ub!=9890:
                    b=np.float(mp.nstr(UInterval[i][0]).split(',')[1][:-1])
                    if b<ub:
                        print("Really?! Max Error!!")
                        model.write("anamoly.lp")
                        print("Interval: ",b)
                        print("Optimization: ",ub)
                U[i][0]=mp.mpi(lb,ub)
        else:
            for i in range(self.n):
                lb=U[i][0][0]
                ub=U[i][0][1]
                #print(lb,ub)
                U[i][0]=mp.mpi(lb,ub)
        #print("Ret>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        return U

    @staticmethod
    def compacify(ORS,Theta,Ac,t):
        n=Ac.shape[0]
        C=ORS[0] # The center is always assumed to be 0 as of now
        V=ORS[1]
        P=ORS[2]

        '''print(C)
        print(V)
        print(P)
        print(t)
        print(Ac)'''

        sv=V.shape[0]
        aS=V.shape[1]
        n=Ac.shape[0]

        U=np.zeros((n,1),dtype=object)

        for i in range(sv):
            s=0
            for j in range(aS):
                s=s+(mp.mpi(P[j][0],P[j][1])*V[i][j])
            s=s+C[i]
            s_min=float(mp.nstr(s).split(',')[0][1:])
            s_max=float(mp.nstr(s).split(',')[1][:-1])
            U[i][0]=(s_min,s_max)

        #print(U)

        Ac=LA.matrix_power(Ac,t)
        C_new=Theta[0] # The center is always assumed to be 0 as of now
        V_new=np.matmul(Ac,Theta[1])
        P_new=[]

        semiDefFlag=False
        modelNewMax = Model("qp")
        modelNewMax.setParam('OutputFlag', False )

        # Create Reachable Set Variables
        predVars=[]
        for i in range(n):
            name="Pred"+str(i)
            #predVars.append(modelNewMax.addVar(-5,5,name=name,vtype='C'))
            predVars.append(modelNewMax.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(n):
            con=0
            for j in range(n):
                con=con+(V_new[i][j]*predVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            #modelNew.addConstr(con<=U[i][0][0],name+".1")
            modelNewMax.addConstr(con==U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(n):
            obj=obj+(predVars[i])
        modelNewMax.setObjective(obj,GRB.MINIMIZE)
        #----------------------------------------

        # Solve the Optimization Problem
        UMax=np.zeros((n,1))
        try:
            modelNewMax.optimize()
            modelNewMax.write('dump.lp')
            status = modelNewMax.Status
            if status==GRB.Status.UNBOUNDED:
                print("UNBOUNDED")
            else:
                if status == GRB.Status.INF_OR_UNBD or \
                   status == GRB.Status.INFEASIBLE  or \
                   status == GRB.Status.UNBOUNDED:
                    print('**The 5model cannot be solved because it is infeasible or unbounded**')
                else:
                    k=0
                    for v in modelNewMax.getVars():
                        #print('%s %g' % (v.varName, v.x))
                        UMax[k][0]=v.x
                        k=k+1

        except:
            semiDefFlag=True

        #print(UMax)
        #exit(0)

        semiDefFlag=False
        modelNewMin = Model("qp")
        modelNewMin.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")


        # Create Reachable Set Variables
        predVars=[]
        for i in range(n):
            name="Pred"+str(i)
            predVars.append(modelNewMin.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(n):
            con=0
            for j in range(n):
                con=con+(V_new[i][j]*predVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            modelNewMin.addConstr(con==U[i][0][0],name+".1")
            #modelNewMin.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(n):
            obj=obj+(predVars[i])
        modelNewMin.setObjective(obj,GRB.MAXIMIZE)
        #----------------------------------------

        # Solve the Optimization Problem
        UMin=np.zeros((n,1))
        try:
            modelNewMin.optimize()
            modelNewMin.write('dump.lp')
            status = modelNewMin.Status
            if status==GRB.Status.UNBOUNDED:
                print("UNBOUNDED")
            else:
                if status == GRB.Status.INF_OR_UNBD or \
                   status == GRB.Status.INFEASIBLE  or \
                   status == GRB.Status.UNBOUNDED:
                    print('**The model cannot be solved because it is infeasible or unbounded**')
                else:
                    k=0
                    for v in modelNewMin.getVars():
                        #print('%s %g' % (v.varName, v.x))
                        UMin[k][0]=v.x
                        k=k+1

        except:
            semiDefFlag=True


        for i in range(n):
            #print((UMin[i][0],UMax[i][0]))
            P_new.append((min(UMin[i][0],UMax[i][0]),max(UMin[i][0],UMax[i][0])))

        #print(P_new)
        #exit(0)

        ORS_new=(C_new,V_new,P_new)

        #print(ORS_new)
        #print("B")
        #exit(0)

        return ORS_new


    @staticmethod
    def randomVec(V):
        n=V.shape[0]
        V_new=np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                V_new[i][j]=random.randint(-4,4)

        return V_new

    @staticmethod
    def filter(V,P):
        n=len(P)
        dim=V.shape[0]
        z=[]
        for i in range(n):
            if abs(P[i][0]-P[i][1])<PRED_EP:
                z.append(i)
        V_new=np.delete(V,z,axis=1)

        sz=V_new.shape[1]
        if sz>dim:
            V_new=np.delete(V_new,list(range(dim,sz)),axis=1)
        #print(V_new.shape)
        #print(V_new)
        V_new=np.identity(dim)
        return V_new

    @staticmethod
    def star2Zono(RS):
        C=RS[0]
        V=RS[1]
        P=RS[2]

        r=V.shape[0]
        c=V.shape[1]

        x=np.asarray(C).reshape(r,1)

        G=np.zeros(V.shape)
        for i in range(r):
            for j in range(c):
                G[i][j]=V[i][j]*P[j][1]

        Z=zonotope(x,G)

        return Z

    @staticmethod
    def zono2Star(Z):
        C=Z.x
        V=Z.G
        P=[(-1,1)]
        n=V.shape[0]
        c=V.shape[1]
        C=C.reshape(1,n)[0]
        P=P*c
        RS=(C,V,P)

        return RS


    @staticmethod
    def lowerDimProj(ORS):

        C=ORS[0]
        V=ORS[1]
        P=ORS[2]

        n=V.shape[0]

        v_len=V.shape[1]

        '''if v_len<n*VLENGTH:
            return ORS'''

        #print("PROJECTING ",n,v_len)

        #BV=Split.filter(V,P)
        Z=Split.star2Zono(ORS)
        G_r=Split.filter(V,P)

        #print("---")
        #print("\n",Z.x)
        #print("\n",Z.G)
        #print(G_r)

        Z_r,Z_r_list,eps=zonotope_order_reduction_outer(Z,G_r,delta=0.01,scale=1.06)
        i_max=len(eps)

        #ORS_compact=Split.zono2Star(Z_r_list[i_max-1])
        ORS_compact=Split.zono2Star(Z_r)

        return ORS_compact

    @staticmethod
    def appxSB(ORS):
        C=ORS[0] # The center is always assumed to be 0 as of now
        V=ORS[1]
        P=ORS[2]

        sv=V.shape[0]
        aS=V.shape[1]

        U=np.zeros((sv,1),dtype=object)

        for i in range(sv):
            s=0
            for j in range(aS):
                s=s+(mp.mpi(P[j][0],P[j][1])*V[i][j])
            s=s+C[i]
            s_min=float(mp.nstr(s).split(',')[0][1:])
            s_max=float(mp.nstr(s).split(',')[1][:-1])
            U[i][0]=(s_min,s_max)

        #print(U)
        C_new=C
        V_new=np.identity(sv)
        P_new=[]

        for i in range(sv):
            P_new.append((U[i][0][0],U[i][0][1]))

        starNew=(C_new,V_new,P_new)

        return starNew

    def getReachableSet(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        #intervalPlot=math.ceil(self.T/7)
        intervalPlot=INTERVAL
        lPlots=[]
        add_stars=0
        prod_mat_stars=0
        comp_U=0
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        ORS_old=self.Theta
        RS=self.Theta
        SRS=[self.Theta]
        U=cu.computeUI_Interval(ORS)
        t=1

        '''(X,Y)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[([(X,Y)],X1,Y1,X2,Y2,X3,Y3)]
        lPlots3=(X2,Y2,X3,Y3)
        #Visualization.displayPlot(s1,s2,lPlots,name+"_0")
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")'''

        while (t<=self.T):
            #sys.stdout.write('\r')
            #sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            #sys.stdout.flush()
            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatStars(SRS)
            ORS_old=ORS
            U_old=U

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)

            '''if t%intervalPlot==0:
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(lst,X1,Y1,X2,Y2,X3,Y3)]
                lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                #Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)'''


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            t=t+1
        #print("\n")
        #time_taken=time.time()-start_time
        #print("Time Taken: ",time_taken)
        #print("")
        return ORS

    def getReachableSetAllList(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        #intervalPlot=math.ceil(self.T/7)
        intervalPlot=INTERVAL
        lPlots=[]
        add_stars=0
        prod_mat_stars=0
        comp_U=0
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        ORS_old=self.Theta
        RS=self.Theta
        SRS=[self.Theta]
        U=cu.computeUI_Interval(ORS)
        t=1

        ORS_List=[ORS]
        RS_List=[RS]

        '''(X,Y)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[([(X,Y)],X1,Y1,X2,Y2,X3,Y3)]
        lPlots3=(X2,Y2,X3,Y3)
        #Visualization.displayPlot(s1,s2,lPlots,name+"_0")
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")'''

        while (t<=self.T):
            #sys.stdout.write('\r')
            #sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            #sys.stdout.flush()
            RS=CompU.prodMatStars(self.A,RS)
            RS_List.append(RS)

            SRS=sample.prodMatStars(SRS)
            ORS_old=ORS
            U_old=U

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            ORS_List.append(ORS)

            '''if t%intervalPlot==0:
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(lst,X1,Y1,X2,Y2,X3,Y3)]
                lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                #Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)'''


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            t=t+1
        #print("\n")
        #time_taken=time.time()-start_time
        #print("Time Taken: ",time_taken)
        #print("")
        return (ORS_List,RS_List)

    def printReachableSetCompactOld(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #ORS_old=self.Theta
        RS=self.Theta
        ORS_compact=self.Theta
        U=cu.computeUI_Interval(ORS)
        U_compact=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
        #lPlots3=(X2,Y2,X3,Y3)
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            RS=CompU.prodMatStars(self.A,RS)

            #SRS=sample.prodMatStars(SRS)
            #ORS_old=ORS
            #U_old=U

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            ORS_compact=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact),U_compact)
            ORS_compact=Split.compacify(ORS_compact,self.Theta,self.Ac,t)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                #(X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
                #lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                #nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                #Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            U_compact=cu.computeUI_Interval(ORS_compact)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time
        print("Time Taken: ",time_taken)
        print("")

    def printReachableSetCompactOld2(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        add_stars=0
        prod_mat_stars=0
        comp_U=0
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #ORS_old=self.Theta
        #RS=self.Theta
        #ORS_compact=self.Theta
        U=cu.computeUI_Interval(ORS)
        #U_compact=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        '''(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
        #lPlots3=(X2,Y2,X3,Y3)
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")'''
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")

        while (t<=self.T):
            '''sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()'''
            #RS=CompU.prodMatStars(self.A,RS)

            #SRS=sample.prodMatStars(SRS)
            #ORS_old=ORS
            #U_old=U

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            #ORS_compact=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact),U_compact)
            #ORS_compact=Split.filterPred(ORS_compact)

            '''if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                #(X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
                #lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                #nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                #Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)'''


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            #U_compact=cu.computeUI_Interval(ORS_compact)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time
        print("Total Time Taken: %.2f" % time_taken)
        print("")
        print("Profiling Results")
        print("- Addition of Stars: %.4f (%.2f %%)" % (Profiling.add_stars, (Profiling.add_stars*100)/time_taken))
        print("- Product of Matrix and Stars: %.4f (%.2f %%)" % (Profiling.prod_mat_stars, (Profiling.prod_mat_stars*100)/time_taken))
        print("- Generating Predicates while computing U (Interval Arithmetic): %.4f (%.2f %%)" % (Profiling.comp_U,(Profiling.comp_U*100)/time_taken))

    def printReachableSetCompactOld3(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        add_stars=0
        prod_mat_stars=0
        comp_U=0
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #ORS_old=self.Theta
        RS=self.Theta
        ORS_compact=self.Theta
        U=cu.computeUI_Interval(ORS)
        U_compact=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        '''(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
        #lPlots3=(X2,Y2,X3,Y3)
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")'''
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)
            #SRS=sample.prodMatStars(SRS)
            #ORS_old=ORS
            #U_old=U

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            ORS_compact=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact),U_compact)
            ORS_compact=Split.lowerDimProj(ORS_compact)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                #(X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
                #lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                #nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                #Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            U_compact=cu.computeUI_Interval(ORS_compact)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time

    def printReachableSetCompactTimeOld(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        time_taken_interval=0
        time_taken_reduction=0
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #ORS_old=self.Theta
        RS=self.Theta
        ORS_compact=self.Theta
        t_reduction=0
        t_interval=time.time()
        U=cu.computeUI_Interval(ORS)
        time_taken_interval=time.time()-t_interval
        time_taken_reduction=time_taken_interval
        U_compact=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        '''(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
        #lPlots3=(X2,Y2,X3,Y3)
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")'''
        #Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            t_interval=time.time()
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)
            time_taken_interval=time_taken_interval+(time.time()-t_interval)

            t_reduction=time.time()
            ORS_compact=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact),U_compact)
            ORS_compact=Split.lowerDimProj(ORS_compact)
            U_compact=cu.computeUI_Interval(ORS_compact)
            time_taken_reduction=time_taken_reduction+(time.time()-t_reduction)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
                name=n+"_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)

            t=t+1

        print("\n")
        print("----Timining Details for "+n+"----")
        print("Without Reduction: ",time_taken_interval)
        print("With Reduction: ",time_taken_reduction)

    def printReachableSetCompactTime(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        time_taken_interval=0
        time_taken_reduction=0
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #ORS_old=self.Theta
        RS=self.Theta
        ORS_compact=self.Theta
        t_reduction=0
        t_interval=time.time()
        U=cu.computeUI_Interval(ORS)
        time_taken_interval=time.time()-t_interval
        time_taken_reduction=time_taken_interval
        U_compact=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            t_interval=time.time()
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)
            time_taken_interval=time_taken_interval+(time.time()-t_interval)

            t_reduction=time.time()
            ORS_compact=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact),U_compact)
            if t%RED_INT==0:
                ORS_compact=Split.appxSB(ORS_compact)
                #ORS_compact=Split.lowerDimProj(ORS_compact)
            U_compact=cu.computeUI_Interval(ORS_compact)
            time_taken_reduction=time_taken_reduction+(time.time()-t_reduction)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lPlots=[(X1,Y1,X2,Y2,X3,Y3)]
                name=n+"_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)

            t=t+1

        print("\n")
        print("----Timining Details for "+n+"----")
        print("Without Reduction: ",time_taken_interval)
        print("With Reduction: ",time_taken_reduction)
        print("Total Time: ",time.time()-start_time_total)

    def printReachableSetMultIS(self,s1,s2,n,Th1,Th2,Th3):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        print()
        print(n)
        print("-----------------\n\n")

        start_time_total=time.time()
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        cu=CompU(self.A,self.Er)
        ORS=self.Theta
        ORS2=Th1
        ORS3=Th2
        ORS4=Th3
        RS=self.Theta
        U=cu.computeUI_Interval(ORS)
        U2=cu.computeUI_Interval(ORS2)
        U3=cu.computeUI_Interval(ORS3)
        U4=cu.computeUI_Interval(ORS4)
        t=1

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS2).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS3).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS4).getPlotsLineFine()
        lPlots=[X1,Y1,X2,Y2,X3,Y3,X4,Y4]
        name=n+"_"+str(t)
        Visualization.displayPlotMultIS(s1,s2,lPlots,name)

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)

            ORS2=CompU.addStars(CompU.prodMatStars(self.Ac,ORS2),U2)
            U2=cu.computeUI_Interval(ORS2)

            ORS3=CompU.addStars(CompU.prodMatStars(self.Ac,ORS3),U3)
            U3=cu.computeUI_Interval(ORS3)

            ORS4=CompU.addStars(CompU.prodMatStars(self.Ac,ORS4),U4)
            U4=cu.computeUI_Interval(ORS4)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS2).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS3).getPlotsLineFine()
                (X4,Y4)=Visualization(s1,s2,ORS4).getPlotsLineFine()
                lPlots=[X1,Y1,X2,Y2,X3,Y3,X4,Y4]
                name=n+"_"+str(t)
                Visualization.displayPlotMultIS(s1,s2,lPlots,name)

            t=t+1

        print("\n")
        print("----Timining Details for "+n+"----")
        print("Total Time: ",time.time()-start_time_total)

    def printReachableSetGrid(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        print()
        print(n)
        print("-----------------\n\n")
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        ORS_o=self.Theta
        ORS_g=Grid(self.Theta,s1,s2).splitStar()
        U_o=cu.computeUI_Interval(ORS_o)
        U_g=[]
        for rs in ORS_g:
            U_g.append(cu.computeUI_Interval(rs))
        t=1

        plotList_g=Visualization.getPlotsLineFineList(ORS_g,s1,s2)
        (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
        Visualization.displayPlotList(s1,s2,plotList_g,(X,Y),name+"_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            ORS_o=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_o),U_o)
            ORS_g=CompU.addStarsList(CompU.prodMatStarsList(self.Ac,ORS_g),U_g)

            if t%intervalPlot==0:
                plotList_g=Visualization.getPlotsLineFineList(ORS_g,s1,s2)
                (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
                name=n+"_"+str(t)
                Visualization.displayPlotList(s1,s2,plotList_g,(X,Y),name)

            U_o=cu.computeUI_Interval(ORS_o)
            tmp=[]
            for rs in ORS_g:
                tmp.append(cu.computeUI_Interval(rs))
            U_g=copy.copy(tmp)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time

    def printReachableSetGridMat(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        images=[]
        intervalPlot=INTERVAL
        lPlots=[]
        lPlotsList=[]
        print()
        print(n)
        print("-----------------\n\n")
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        ORS_o=self.Theta
        ORS_g1=self.Theta
        ORS_g2=self.Theta
        ORS_g3=self.Theta
        ORS_g4=self.Theta
        erList=GridMat(self.A,self.Er).splitEr()
        cu1=CompU(self.A,erList[0])
        cu2=CompU(self.A,erList[1])
        cu3=CompU(self.A,erList[2])
        cu4=CompU(self.A,erList[3])
        #print(erList[2],self.Er)
        Ac1=cu1.computeCenter()
        Ac2=cu2.computeCenter()
        Ac3=cu3.computeCenter()
        Ac4=cu4.computeCenter()
        U_o=cu.computeUI_Interval(ORS_o)
        U_g1=cu1.computeUI_Interval(ORS_g1)
        U_g2=cu2.computeUI_Interval(ORS_g2)
        U_g3=cu3.computeUI_Interval(ORS_g3)
        U_g4=cu1.computeUI_Interval(ORS_g4)

        print("----Eigenvalues----")
        print("Ac: ",LA.eig(self.Ac)[0])
        print(LA.eig(Ac1)[0])
        print(LA.eig(Ac2)[0])
        print(LA.eig(Ac3)[0])
        print(LA.eig(Ac4)[0])
        print("------\n\n")
        print("----Singular Values----")
        print("Ac: ",LA.norm(self.Ac,2))
        print(LA.norm(Ac1,2))
        print(LA.norm(Ac2,2))
        print(LA.norm(Ac3,2))
        print(LA.norm(Ac4,2))
        print("------\n\n")
        t=1
        (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
        (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS_g3).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS_g4).getPlotsLineFine()
        Visualization.displayPlotList(s1,s2,[(X1,Y1),(X2,Y2),(X3,Y3),(X4,Y4)],(X,Y),name+"_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            ORS_o=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_o),U_o)
            ORS_g1=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g1),U_g1)
            ORS_g2=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g2),U_g2)
            ORS_g3=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g3),U_g3)
            ORS_g4=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g4),U_g4)

            if t%intervalPlot==0:
                (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
                (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS_g3).getPlotsLineFine()
                (X4,Y4)=Visualization(s1,s2,ORS_g4).getPlotsLineFine()

                name=n+"_"+str(t)
                images.append(Visualization.getPlotList(s1,s2,[(X1,Y1),(X2,Y2),(X3,Y3),(X4,Y4)],(X,Y),name))

            U_o=cu.computeUI_Interval(ORS_o)
            U_g1=cu.computeUI_Interval(ORS_g1)
            U_g2=cu.computeUI_Interval(ORS_g2)
            U_g3=cu.computeUI_Interval(ORS_g3)
            U_g4=cu.computeUI_Interval(ORS_g4)
            t=t+1
        print("\n")
        images[0].save("TempGIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=30, loop=0)
        time_taken=time.time()-start_time

    def printReachableSetGridOneMat(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        images=[]
        intervalPlot=INTERVAL
        lPlots=[]
        lPlotsList=[]
        print()
        print(n)
        print("-----------------\n\n")
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        ORS_o=self.Theta
        ORS_g1=self.Theta
        ORS_g2=self.Theta
        erList=GridMat(self.A,self.Er).splitErOne()
        cu1=CompU(self.A,erList[0])
        cu2=CompU(self.A,erList[1])
        Ac1=cu1.computeCenter()
        Ac2=cu2.computeCenter()
        U_o=cu.computeUI_Interval(ORS_o)
        U_g1=cu1.computeUI_Interval(ORS_g1)
        U_g2=cu2.computeUI_Interval(ORS_g2)

        print("----Eigenvalues----")
        print("Ac: ",LA.eig(self.Ac)[0])
        print(LA.eig(Ac1)[0])
        print(LA.eig(Ac2)[0])
        print("------\n\n")
        print("----Singular Values----")
        print("Ac: ",LA.norm(self.Ac,2))
        print(LA.norm(Ac1,2))
        print(LA.norm(Ac2,2))
        print("------\n\n")
        t=1
        (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
        (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()
        Visualization.displayPlotList(s1,s2,[(X1,Y1),(X2,Y2),([],[]),([],[])],(X,Y),name+"_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            ORS_o=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_o),U_o)
            ORS_g1=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g1),U_g1)
            ORS_g2=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_g2),U_g2)

            if t%intervalPlot==0:
                (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
                (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()

                name=n+"_"+str(t)
                images.append(Visualization.getPlotList(s1,s2,[(X1,Y1),(X2,Y2),([],[]),([],[])],(X,Y),name))

            U_o=cu.computeUI_Interval(ORS_o)
            U_g1=cu.computeUI_Interval(ORS_g1)
            U_g2=cu.computeUI_Interval(ORS_g2)
            t=t+1
        print("\n")
        images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        time_taken=time.time()-start_time

    def printReachableSetGridOneMat2(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        images=[]
        intervalPlot=INTERVAL
        lPlots=[]
        lPlotsList=[]
        print()
        print(n)
        print("-----------------\n\n")
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        ORS_o=self.Theta
        ORS_g1=self.Theta
        ORS_g2=self.Theta
        ORS_g3=self.Theta
        ORS_g4=self.Theta
        erList=GridMat(self.A,self.Er).splitErOne2()
        cu1=CompU(self.A,erList[0])
        cu2=CompU(self.A,erList[1])
        cu3=CompU(self.A,erList[2])
        cu4=CompU(self.A,erList[3])
        #print(erList[2],self.Er)
        Ac1=cu1.computeCenter()
        Ac2=cu2.computeCenter()
        Ac3=cu3.computeCenter()
        Ac4=cu4.computeCenter()
        U_o=cu.computeUI_Interval(ORS_o)
        U_g1=cu1.computeUI_Interval(ORS_g1)
        U_g2=cu2.computeUI_Interval(ORS_g2)
        U_g3=cu3.computeUI_Interval(ORS_g3)
        U_g4=cu1.computeUI_Interval(ORS_g4)

        print("----Eigenvalues----")
        print("Ac: ",LA.eig(self.Ac)[0])
        print(LA.eig(Ac1)[0])
        print(LA.eig(Ac2)[0])
        print(LA.eig(Ac3)[0])
        print(LA.eig(Ac4)[0])
        print("------\n\n")
        print("----Singular Values----")
        print("Ac: ",LA.norm(self.Ac,2))
        print(LA.norm(Ac1,2))
        print(LA.norm(Ac2,2))
        print(LA.norm(Ac3,2))
        print(LA.norm(Ac4,2))
        print("------\n\n")
        t=1
        (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
        (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS_g3).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS_g4).getPlotsLineFine()
        Visualization.displayPlotList(s1,s2,[(X1,Y1),(X2,Y2),(X3,Y3),(X4,Y4)],(X,Y),name+"_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            ORS_o=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_o),U_o)
            ORS_g1=CompU.addStars(CompU.prodMatStars(Ac1,ORS_g1),U_g1)
            ORS_g2=CompU.addStars(CompU.prodMatStars(Ac2,ORS_g2),U_g2)
            ORS_g3=CompU.addStars(CompU.prodMatStars(Ac3,ORS_g3),U_g3)
            ORS_g4=CompU.addStars(CompU.prodMatStars(Ac4,ORS_g4),U_g4)

            if t%intervalPlot==0:
                (X,Y)=Visualization(s1,s2,ORS_o).getPlotsLineFine()
                (X1,Y1)=Visualization(s1,s2,ORS_g1).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_g2).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS_g3).getPlotsLineFine()
                (X4,Y4)=Visualization(s1,s2,ORS_g4).getPlotsLineFine()

                name=n+"_"+str(t)
                images.append(Visualization.getPlotList(s1,s2,[(X1,Y1),(X2,Y2),(X3,Y3),(X4,Y4)],(X,Y),name))

            U_o=cu.computeUI_Interval(ORS_o)
            U_g1=cu.computeUI_Interval(ORS_g1)
            U_g2=cu.computeUI_Interval(ORS_g2)
            U_g3=cu.computeUI_Interval(ORS_g3)
            U_g4=cu.computeUI_Interval(ORS_g4)
            t=t+1
        print("\n")
        images[0].save("TempGIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=30, loop=0)
        time_taken=time.time()-start_time

    def printReachableSetAll(self,s1,s2,flow,n):
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        images=[]
        time_taken_interval=0
        time_taken_reduction=0
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        SRS=[self.Theta]
        RS=self.Theta
        ORS_compact_zono=self.Theta
        ORS_compact_interval=self.Theta
        t_reduction_zono=0
        t_reduction_interval=0
        t_interval=time.time()
        U=cu.computeUI_Interval(ORS)
        time_taken_interval=time.time()-t_interval
        time_taken_reduction_zono=time_taken_interval
        time_taken_reduction_interval=time_taken_interval
        U_compact_zono=U
        U_compact_interval=U
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        #(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        #(X2,Y2)=Visualization(s1,s2,ORS_compact_zono).getPlotsLineFine()
        #(X3,Y3)=Visualization(s1,s2,ORS_compact_interval).getPlotsLineFine()
        #(X4,Y4)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #(X5,Y5)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        #lPlots=[[(X5,Y5)],[X1,Y1,X2,Y2,X3,Y3,X4,Y4]]
        #images.append(Visualization.getPlotAll(s1,s2,lPlots,name+"_0"))

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatStars(SRS)

            t_interval=time.time()
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)
            time_taken_interval=time_taken_interval+(time.time()-t_interval)

            t_reduction_interval=time.time()
            ORS_compact_interval=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact_interval),U_compact_interval)
            if t%RED_INT_INTRVL==0:
                ORS_compact_interval=Split.appxSB(ORS_compact_interval)
            U_compact_interval=cu.computeUI_Interval(ORS_compact_interval)
            time_taken_reduction_interval=time_taken_reduction_interval+(time.time()-t_reduction_interval)

            t_reduction_zono=time.time()
            ORS_compact_zono=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact_zono),U_compact_zono)
            if t%RED_INT_ZONO==0:
                ORS_compact_zono=Split.lowerDimProj(ORS_compact_zono)
            U_compact_zono=cu.computeUI_Interval(ORS_compact_zono)
            time_taken_reduction_zono=time_taken_reduction_zono+(time.time()-t_reduction_zono)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact_zono).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS_compact_interval).getPlotsLineFine()
                (X4,Y4)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                lPlots=[lst,[X1,Y1,X2,Y2,X3,Y3,X4,Y4]]
                name=n+"_"+str(t)
                images.append(Visualization.getPlotAll(s1,s2,lPlots,flow,name))

            t=t+1

        print("\n")
        #images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        print("----Timining Details for "+n+"----")
        print("Without Reduction: ",time_taken_interval)
        print("With Reduction (Zono): ",time_taken_reduction_zono)
        print("With Reduction (Interval): ",time_taken_reduction_interval)
        print("Total Time: ",time.time()-start_time_total)

    def getReachableSetAll(self,s1=0,s2=0,n="RS"):
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        images=[]
        time_taken_interval=0
        time_taken_reduction=0
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #SRS=[self.Theta]
        #RS=self.Theta
        #ORS_compact_zono=self.Theta
        #ORS_compact_interval=self.Theta
        #t_reduction_zono=0
        #t_reduction_interval=0
        #t_interval=time.time()
        U=cu.computeUI_Interval(ORS)
        #time_taken_interval=time.time()-t_interval
        #time_taken_reduction_zono=time_taken_interval
        #time_taken_reduction_interval=time_taken_interval
        #U_compact_zono=U
        #U_compact_interval=U
        t=1
        print()
        #print(n)
        print("-----------------\n\n")

        #(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        #(X2,Y2)=Visualization(s1,s2,ORS_compact_zono).getPlotsLineFine()
        #(X3,Y3)=Visualization(s1,s2,ORS_compact_interval).getPlotsLineFine()
        #(X4,Y4)=Visualization(s1,s2,ORS).getPlotsLineFine()
        #(X5,Y5)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        #lPlots=[[(X5,Y5)],[X1,Y1,X2,Y2,X3,Y3,X4,Y4]]
        #images.append(Visualization.getPlotAll(s1,s2,lPlots,name+"_0"))

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            #RS=CompU.prodMatStars(self.A,RS)

            #SRS=sample.prodMatStars(SRS)

            t_interval=time.time()
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)
            time_taken_interval=time_taken_interval+(time.time()-t_interval)

            '''t_reduction_interval=time.time()
            ORS_compact_interval=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact_interval),U_compact_interval)
            if t%RED_INT_INTRVL==0:
                ORS_compact_interval=Split.appxSB(ORS_compact_interval)
            U_compact_interval=cu.computeUI_Interval(ORS_compact_interval)
            time_taken_reduction_interval=time_taken_reduction_interval+(time.time()-t_reduction_interval)
            '''
            '''
            t_reduction_zono=time.time()
            ORS_compact_zono=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_compact_zono),U_compact_zono)
            if t%RED_INT_ZONO==0:
                ORS_compact_zono=Split.lowerDimProj(ORS_compact_zono)
            U_compact_zono=cu.computeUI_Interval(ORS_compact_zono)
            time_taken_reduction_zono=time_taken_reduction_zono+(time.time()-t_reduction_zono)
            '''

            if False:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS_compact_zono).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS_compact_interval).getPlotsLineFine()
                (X4,Y4)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                lPlots=[lst,[X1,Y1,X2,Y2,X3,Y3,X4,Y4]]
                name=n+"_"+str(t)
                images.append(Visualization.getPlotAll(s1,s2,lPlots,name))

            t=t+1

        print("\n")
        t=t-1
        #images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        print("----Timining Details for "+n+"----")
        print("Without Reduction: ",time_taken_interval)
        #print("With Reduction (Zono): ",time_taken_reduction_zono)
        #print("With Reduction (Interval): ",time_taken_reduction_interval)
        #print("Total Time: ",time.time()-start_time_total)

        '''(X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS_compact_zono).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS_compact_interval).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS).getPlotsLineFine()
        lst=Sampling.getPlotsLineFine(s1,s2,SRS)
        lPlots=[lst,[X1,Y1,X2,Y2,X3,Y3,X4,Y4]]
        name=n+"_"+str(t)
        images.append(Visualization.getPlotAll(s1,s2,lPlots,name))'''

        return ORS
        #return (ORS,ORS_compact_zono,ORS_compact_interval)

    def printReachableSetPred(self,s1,s2,pnew,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        nameU=n
        #intervalPlot=math.ceil(self.T/7)
        intervalPlot=INTERVAL
        lPlots=[]
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        ORS_old=self.Theta
        RS=self.Theta
        SRS=[self.Theta]
        U=cu.computeUI_Pred(ORS,pnew)
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X,Y)=Visualization(s1,s2,SRS[0]).getPlotsLineFinePred()
        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFinePred()
        (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFinePred()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFinePred()
        #lPlots.append((X1,Y1,X2,Y2,X3,Y3))
        lPlots=[([(X,Y)],X1,Y1,X2,Y2,X3,Y3)]
        lPlots3=(X2,Y2,X3,Y3)
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")
        Visualization.displayPlotSingle(s1,s2,lPlots3,nameU+"U_0")

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatStars(SRS)
            ORS_old=ORS
            U_old=U

            ORS=CompU.addStarsPred(CompU.prodMatStars(self.Ac,ORS),U)
            #print(ORS)
            #exit(0)

            if t%intervalPlot==0:
                lst=Sampling.getPlotsLineFinePred(s1,s2,SRS)
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFinePred()
                (X2,Y2)=Visualization(s1,s2,ORS_old).getPlotsLineFinePred()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFinePred()
                lPlots=[(lst,X1,Y1,X2,Y2,X3,Y3)]
                lPlots2=(X2,Y2,X3,Y3)
                name=n+"_"+str(t)
                nameU=n+"U_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)
                Visualization.displayPlotSingle(s1,s2,lPlots2,nameU)


            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Pred(ORS,pnew)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time
        print("Time Taken: ",time_taken)
        print("")

    @staticmethod
    def getVectorLength(RS):
        return RS[1].shape[1]

    def reachableSetProfile(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''

        print()
        print(n)
        print("-----------------\n\n")

        name=n
        nameU=n
        strRep=""
        intervalPlot=10
        lPlots=[]
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        ORS=self.Theta
        RS=self.Theta

        U=cu.computeUI_Interval(ORS)
        t=1

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            RS=CompU.prodMatStars(self.A,RS)
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            start_time_interval=time.time()
            U=cu.computeUI_Interval(ORS)
            time_taken_interval=time.time()-start_time_interval
            if t%intervalPlot==0:
                ORS_compact=ORS
                if Split.getVectorLength(ORS)<900:
                    start_time_reduction=time.time()
                    ORS_compact=Split.lowerDimProj(ORS)
                    time_taken_reduction=time.time()-start_time_reduction
                    (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                    (X2,Y2)=Visualization(s1,s2,ORS_compact).getPlotsLineFine()
                    lPlots=[X1,Y1,X2,Y2]
                    name=n+"_"+str(t)
                    Visualization.displayPlotSingle(s1,s2,lPlots,name)
                else:
                    time_taken_reduction=0
                #print("Vector Size: ",Split.getVectorLength(ORS),"; Interval: ",time_taken_interval,"; Reduction: ",time_taken_reduction)
                strRep=strRep+"Vector Size: "+str(Split.getVectorLength(ORS))+"; Interval: "+str(time_taken_interval)+"; Reduction: "+str(time_taken_reduction)+"\n"
            t=t+1
        print("\n")
        time_taken=time.time()-start_time
        print("\n\n----------Profiling Report for "+n+"----------")
        print(strRep)
        print("\nTotal Time: ",time_taken)

    def printStarPredReport(RS,i):
        pred=RS[2]
        n=len(pred)
        z=0

        for p in pred:
            if abs(p[0]-p[1])<EPSILON:
                z=z+1
        print("At Interval ",i,": ",z,"/",n," (",((z*100)/n),"% )")

    def filterPred(RS):
        C=RS[0]
        V=RS[1]
        P=RS[2]
        z=[]
        P_new=[]
        n=len(P)
        for i in range(n):
            if abs(P[i][0]-P[i][1])<EPSILON:
                z.append(i)
            else:
                P_new.append(P[i])
        V_new=np.delete(V,z,axis=1)

        return (C,V_new,P_new)

    def printReachableSetRand(self,s1,s2,n):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        name=n
        intervalPlot=INTERVAL
        lPlots=[]
        print()
        print(n)
        print("-----------------\n\n")
        start_time=time.time()
        cu=CompU(self.A,self.Er)
        #sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        #RS=self.Theta
        #SRS=[self.Theta]
        #ORS_rand=[self.Theta]*SAMPLES
        U=cu.computeUI_Interval(ORS)
        t=1
        #U_rand=cu.computeUI_IntervalRand(ORS_rand)

        '''(X,Y)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS_rand[0]).getPlotsLineFine()
        lPlots=[([(X,Y)],X2,Y2,X3,Y3,[(X4,Y4)])]
        Visualization.displayPlot(s1,s2,lPlots,name+"_0")'''

        while (t<=self.T):
            '''sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()'''
            #RS=CompU.prodMatStars(self.A,RS)
            #SRS=sample.prodMatStars(SRS)
            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            #ORS_rand=CompU.addStarsList(CompU.prodMatStarsList(self.Ac,ORS_rand),U_rand)

            if t%intervalPlot==0:
                '''lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                (X2,Y2)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS).getPlotsLineFine()
                lst2=Sampling.getPlotsLineFine(s1,s2,ORS_rand)
                lPlots=[(lst,X2,Y2,X3,Y3,lst2)]
                name=n+"_"+str(t)
                Visualization.displayPlot(s1,s2,lPlots,name)'''
                #Split.printStarPredReport(ORS,t)



            '''print("Center of ORS: ",ORS[0])
            print("Vector of ORS: ")
            print(ORS[1])
            print("Predicate of ORS: ",ORS[2])
            print("\n\n")'''


            U=cu.computeUI_Interval(ORS)
            #U_rand=cu.computeUI_IntervalRand(ORS_rand)
            t=t+1
        print("\n")
        time_taken=time.time()-start_time
        print("Time Taken: ",time_taken)
        print("")

    def printReachableSetTwo(self,s1,s2,Er2,n):
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        images=[]
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        cu2=CompU(self.A,Er2)
        sample=Sampling(self.A,self.Er)
        ORS=self.Theta
        ORS2=self.Theta
        RS=self.Theta
        SRS=[self.Theta]
        U=cu.computeUI_Interval(ORS)
        U2=cu.computeUI_Interval(ORS2)
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,ORS2).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        lPlots=[[(X4,Y4)],[X1,Y1,X2,Y2,X3,Y3]]
        images.append(Visualization.getPlotTwo(s1,s2,lPlots,name+"_0"))

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatStars(SRS)

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)

            ORS2=CompU.addStars(CompU.prodMatStars(self.Ac,ORS2),U2)
            U2=cu2.computeUI_Interval(ORS2)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS2).getPlotsLineFine()
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                lPlots=[lst,[X1,Y1,X2,Y2,X3,Y3]]
                name=n+"_"+str(t)
                images.append(Visualization.getPlotTwo(s1,s2,lPlots,name))

            t=t+1

        print("\n")
        images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        print("----Timining Details for "+n+"----")
        print("Total Time: ",time.time()-start_time_total)

    def printReachableSetOrdComp(self,s1,s2,Er2,n):
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        images=[]
        start_time_total=time.time()
        sample=Sampling(self.A,self.Er)
        sample2=Sampling(self.A,Er2)
        RS=self.Theta
        SRS=[self.Theta]
        SRS2=[self.Theta]
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,SRS[0]).getPlotsLineFine()
        (X3,Y3)=Visualization(s1,s2,SRS2[0]).getPlotsLineFine()
        lPlots=[X1,Y1,[(X2,Y2)],[(X3,Y3)]]
        images.append(Visualization.getPlotOrdComp(s1,s2,lPlots,name+"_0"))

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatStars(SRS)
            SRS2=sample2.prodMatStars(SRS2)

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                lst=Sampling.getPlotsLineFine(s1,s2,SRS)
                lst2=Sampling.getPlotsLineFine(s1,s2,SRS2)
                lPlots=[X1,Y1,lst,lst2]
                name=n+"_"+str(t)
                images.append(Visualization.getPlotOrdComp(s1,s2,lPlots,name))

            t=t+1

        print("\n")
        images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        print("----Timining Details for "+n+"----")
        print("Total Time: ",time.time()-start_time_total)

    def egHeuristics():
        '''
        Just to generate examples for the heuritics
        '''
        C=[0,0]
        V=np.array([
        [0,1,0.5,0.9],
        [1,-1,0.2,0.3]
        ])
        P=[(-10,10),(-10,10),(-2,2),(-1,1)]
        rs=(C,V,P)

        rsInterval=Split.appxSB(rs)
        rsZono=Split.lowerDimProj(rs)

        print(rsInterval)
        print(rsZono)
        Visualization.displayPlotEgHeuristics(0,1,rs,rsInterval,rsZono)


if False:
    A2=np.array([
    [1,1,-2],
    [2,0.2,0],
    [0,0.1,0.1]
    ])
    A=np.array([
    [2,1,-2,1],
    [2,1.2,0,0],
    [0,0.1,1.1,1],
    [0,0,0,1]
    ])
    E={
    (0,1):[0.9,1.1],
    (1,0):[0.8,1.2],
    (2,3):[0.9,1.1]
    }
    E2={
    (0,1):[0.9,1.5],
    (1,0):[0.8,1.2],
    }
    T=200
    IS2=np.array([
    [1],
    [1],
    [1]
    ])
    IS=np.array([
    [1],
    [1],
    [1],
    [1]
    ])
    Split(A,E,IS,T).printReachableSet()

if False:
    A2=np.array([
    [2,1],
    [2,1.2]
    ])
    A=np.array([
    [5,2],
    [2,5]
    ])
    E={
    (0,1): [0.9,1.1]
    }

    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=[(-1,1),(-1,1)]
    rs=(C,V,P)

    T=100

    sp=Split(A,E,rs,T)
    sp.printReachableSetCompact(0,1,"Test")

if False:
    A=np.array([
    [2,1],
    [-2,1.2]
    ])
    E={
    (0,1): [0.98,1.02]
    }

    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=np.array([
    [1,0],
    [0,1]
    ])
    D=np.array([
    [(1,1)],
    [(1,1)]
    ])
    rs=(C,V,(P,D))

    pnew=np.array([
    [1,0],
    [0,1]
    ])

    T=10

    sp=Split(A,E,rs,T)
    sp.printReachableSet(0,1,pnew,"Test")

if False:
    Split.egHeuristics()
