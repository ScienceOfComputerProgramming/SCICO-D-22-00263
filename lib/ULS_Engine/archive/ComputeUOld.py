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

BIGM=1e500
EPSILON=1e-3


class CompU:
    '''
    Computes the reachable set of a given linear discrete dynamical system A
    form the initial set Theta, upto time T
    '''

    def __init__(self,A,E):
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
        self.n=A.shape[0]
        self.Ac=self.computeCenter()
        self.methodName="Compute U"

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

    def computeUOld(self,rs):
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
                        mx=obj.getValue()
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

    def computeU2(self,rs):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2]

        semiDefFlag=False

        model = Model("qp")
        model.setParam('OutputFlag', False )
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
            a=P[i][0]
            b=P[i][1]
            model.addConstr(reachVars[i]<=max(a,b),name+".1")
            model.addConstr(reachVars[i]>=min(a,b),name+".1")
        #---------------------------------

        # Prepare the objective functions and get min max
        U=np.zeros((self.n,1),dtype=object)

        for i in range(self.n):
            obj=0
            for j in range(self.n):
                obj2=0
                for k in range(self.n):
                    if (i,k) in self.Er:
                        pertV=model.getVarByName("Pert"+str((i,k)))
                        obj2=obj2+(((self.A[i][k]*pertV)-self.Ac[i][k])*V[k][j])
                    else:
                        obj2=obj2+((self.A[i][k]-self.Ac[i][k])*V[k][j])
                obj=obj2*reachVars[j]
            objC=0
            for j in range(self.n):
                if (i,j) in self.Er:
                    pertV=model.getVarByName("Pert"+str((i,j)))
                    objC=objC+(((self.A[i][j]*pertV)-self.Ac[i][j])*C[j])
                else:
                    objC=objC+((self.A[i][j]-self.Ac[i][j])*C[j])
            obj=objC+obj

            # Obtain Minimum
            mn=-9890
            model.setObjective(obj,GRB.MINIMIZE)
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
                        mn=obj.getValue()
            except:
                semiDefFlag=True


            mx=9890
            model.setObjective(obj,GRB.MAXIMIZE)
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
                        mx=obj.getValue()
            except:
                semiDefFlag=True


            U[i][0]=(mn,mx)
        #----------------------------------
        print(U)


        # Find a new Star set that overapproximates rs
        semiDefFlag=False
        modelNew = Model("qp")
        modelNew.setParam('OutputFlag', False )

        print("Construction of the new overapproximated star is under construction!")
        C_new=C
        V_new=np.matmul(self.Ac,V)
        P_new=[]

        # Introducing new variables to bloat the Predicate
        predVars=[]
        for i in range(self.n):
            name="Pred"+str(i)
            predVars.append(modelNew.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create Reachable Set Variables
        reachVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            reachVars.append(modelNew.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions
        for i in range(self.n):
            con=0
            for j in range(self.n):
                pertV=model.getVarByName("Pert"+str((i,j)))
                con=con+(V_new[i][j]*reachVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            modelNew.addConstr(con<=U[i][0][0],name+".1")
            modelNew.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Bloat Predicate Constraints
        for i in range(self.n):
            modelNew.addConstr(reachVars[i]<=predVars[i])
            modelNew.addConstr(reachVars[i]>=-predVars[i])
            modelNew.addConstr(predVars[i]==100)
        #---------------------------------------

        # Create the Objective Function
        obj=0
        for i in range(self.n):
            obj=obj+(predVars[i])
        modelNew.setObjective(obj,GRB.MINIMIZE)
        #----------------------------------------

        # Solve the Optimization Problem
        try:
            modelNew.optimize()
            modelNew.write('dump.lp')
            status = modelNew.Status
            if status==GRB.Status.UNBOUNDED:
                print("UNBOUNDED")
            else:
                if status == GRB.Status.INF_OR_UNBD or \
                   status == GRB.Status.INFEASIBLE  or \
                   status == GRB.Status.UNBOUNDED:
                    print('**The model cannot be solved because it is infeasible or unbounded**')
                else:
                    for v in modelNew.getVars():
                        print('%s %g' % (v.varName, v.x))
        except:
            semiDefFlag=True
        #---------------------------------------

        #---------------------------------------------

    def computeU(self,rs):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2]

        '''print(self.A)
        print(self.Ac)
        print(V)'''

        semiDefFlag=False

        model = Model("qp")
        model.setParam('OutputFlag', False )
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
            a=P[i][0]
            b=P[i][1]
            model.addConstr(reachVars[i]<=max(a,b),name+".1")
            model.addConstr(reachVars[i]>=min(a,b),name+".1")
        #---------------------------------

        # Prepare the objective functions and get min max
        U=np.zeros((self.n,1),dtype=object)

        for i in range(self.n):
            obj=0
            for j in range(self.n):
                obj2=0
                for k in range(self.n):
                    if (i,k) in self.Er:
                        pertV=model.getVarByName("Pert"+str((i,k)))
                        obj2=obj2+(((self.A[i][k]*pertV)-self.Ac[i][k])*V[k][j])
                    else:
                        obj2=obj2+((self.A[i][k]-self.Ac[i][k])*V[k][j])
                obj=obj+(obj2*reachVars[j])
            objC=0
            for j in range(self.n):
                if (i,j) in self.Er:
                    pertV=model.getVarByName("Pert"+str((i,j)))
                    objC=objC+(((self.A[i][j]*pertV)-self.Ac[i][j])*C[j])
                else:
                    objC=objC+((self.A[i][j]-self.Ac[i][j])*C[j])
            obj=objC+obj
            #print(obj)

            # Obtain Minimum
            mn=-9890
            model.setObjective(obj,GRB.MINIMIZE)
            try:
                model.optimize()
                #print(obj)
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
                model.write("dumpN.lp")
                print("i: ",i)



            mx=9890
            model.setObjective(obj,GRB.MAXIMIZE)
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
                        mx=obj.getValue()
            except:
                semiDefFlag=True


            U[i][0]=(mn,mx)
        #----------------------------------
        #print(U)
        #exit(0)

        # Find a new Star set that overapproximates rs
        semiDefFlag=False
        modelNewMax = Model("qp")
        modelNewMax.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")
        C_new=C
        V_new=np.matmul(self.Ac,V)
        P_new=[]


        # Create Reachable Set Variables
        reachVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            reachVars.append(modelNewMax.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(self.n):
            con=0
            for j in range(self.n):
                con=con+(V_new[i][j]*reachVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            #modelNew.addConstr(con<=U[i][0][0],name+".1")
            modelNewMax.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(self.n):
            obj=obj+(reachVars[i])
        modelNewMax.setObjective(obj,GRB.MINIMIZE)
        #----------------------------------------

        # Solve the Optimization Problem
        UMax=np.zeros((self.n,1))
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


        semiDefFlag=False
        modelNewMin = Model("qp")
        modelNewMin.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")


        # Create Reachable Set Variables
        reachVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            reachVars.append(modelNewMin.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(self.n):
            con=0
            for j in range(self.n):
                con=con+(V_new[i][j]*reachVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            modelNewMin.addConstr(con<=U[i][0][0],name+".1")
            #modelNewMin.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(self.n):
            obj=obj+(reachVars[i])
        modelNewMin.setObjective(obj,GRB.MAXIMIZE)
        #----------------------------------------

        # Solve the Optimization Problem
        UMin=np.zeros((self.n,1))
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


        for i in range(self.n):
            #print((UMin[i][0],UMax[i][0]))
            P_new.append((UMin[i][0],UMax[i][0]))

        #print(P_new)

        starNew=(C_new,V_new,P_new)

        print("-------Given Over-approximated Star-------")
        print("Center: ",C)
        print("Vector: ")
        print(V)
        print("Predicate (Box): ",P)
        print("-----------------------")
        print()
        print("-------Computed U Star-------")
        print("Center: ",C_new)
        print("Vector: ")
        print(V_new)
        print("Predicate (Box): ",P_new)
        print("-----------------------")
        return starNew
        #-------------------------------------






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
    P=[(5,5),(5,5)]
    rs=(C,V,P)

    #u=CompU(A,E)
    #u.computeU(rs)

    A2=np.array([
    [2,1,-1],
    [-2,1.2,-1.2],
    [-1,-1,1.2]
    ])
    E2={
    (0,1): [0.98,1.02],
    (2,2): [0.98,1.02]
    }

    C2=[0,0,0]
    V2=np.array([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
    P2=[(2,2),(2,2),(2,2)]
    rs2=(C2,V2,P2)

    u2=CompU(A2,E2)
    u2.computeU(rs2)
