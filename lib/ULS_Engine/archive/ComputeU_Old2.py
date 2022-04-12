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
from operator import add

from PredicateGen import *
import Profiling

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

    def computeUI_Pred(self,rs,pnew):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2][0]
        D=rs[2][1]

        sv=V.shape[0]
        aS=V.shape[1]
        n_constr=P.shape[0] # Number of Constraints on the predicate Variables
        #print(P)

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

        # Create Predicate Variables
        predVars=[]
        for i in range(aS):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
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

        # Add Predicate Constraints
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
        #---------------------------------

        # Prepare the objective functions and get min max
        U=np.zeros((self.n,1),dtype=object)

        for i in range(sv):
            obj=0
            for j in range(aS):
                obj2=0
                for k in range(sv):
                    if (i,k) in self.Er:
                        pertV=model.getVarByName("Pert"+str((i,k)))
                        obj2=obj2+(((self.A[i][k]*pertV)-self.Ac[i][k])*V[k][j])
                    else:
                        obj2=obj2+((self.A[i][k]-self.Ac[i][k])*V[k][j])
                obj=obj+(obj2*predVars[j])
            objC=0
            for j in range(sv):
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
                #model.write("dumpN.lp")
                print("Shoot. i: ",i)



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
                print("Shoot. i: ",i)


            U[i][0]=(mn,mx)
        #----------------------------------

        #print(U); exit(0)

        # Find a new Star set that overapproximates rs
        semiDefFlag=False
        modelNewMax = Model("qp")
        modelNewMax.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")
        C_new=np.zeros(self.n)
        V_new=np.identity(self.n)
        P_new=pnew
        D_new=np.zeros((self.n,1),dtype=object)

        # Create Predicate Variables
        predVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            predVars.append(modelNewMax.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(self.n):
            con=0
            for j in range(self.n):
                con=con+(V_new[i][j]*predVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            #modelNew.addConstr(con<=U[i][0][0],name+".1")
            modelNewMax.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(P_new.shape[0]):
            objT=0
            for j in range(self.n):
                objT=objT+(P_new[i][j]*predVars[j])
            obj=obj+(objT)
        modelNewMax.setObjective(obj,GRB.MINIMIZE)
        modelNewMax.optimize()
        #----------------------------------------

        # Solve the Optimization Problem
        UMax=np.zeros((self.n,1))
        try:
            modelNewMax.optimize()
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
            print("Shoot (Finding)")


        semiDefFlag=False
        modelNewMin = Model("qp")
        modelNewMin.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")


        # Create Reachable Set Variables
        predVars=[]
        for i in range(self.n):
            name="IS"+str(i)
            predVars.append(modelNewMin.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Create the star set Constraints functions for Maximum Bound
        for i in range(self.n):
            con=0
            for j in range(self.n):
                con=con+(V_new[i][j]*predVars[j])
            con=C_new[i]+con
            name="Predicate-C-"+str(i)
            modelNewMin.addConstr(con<=U[i][0][0],name+".1")
            #modelNewMin.addConstr(con>=U[i][0][1],name+".2")
        #-------------------------------

        # Create the Objective Function
        obj=0
        for i in range(P_new.shape[0]):
            objT=0
            for j in range(self.n):
                objT=objT+(P_new[i][j]*predVars[j])
            obj=obj+(objT)
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
            print("Shoot (Finding)")


        for i in range(self.n):
            #print((UMin[i][0],UMax[i][0]))
            D_new[i][0]=(UMin[i][0],UMax[i][0])

        #print(P_new)

        starNew=(C_new,V_new,(P_new,D_new))

        '''print("-------Given Over-approximated Star-------")
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
        print("-----------------------")'''
        return starNew
        #-------------------------------------

    def computeUI_IntervalRand(self,rs_list):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        V_randList=PredGen(self.A,self.n,len(rs_list)).getRandPred()

        '''print(V_randList)
        exit(0)'''
        star_list=[]

        for e in range(len(rs_list)):

            rs=rs_list[e]

            C=rs[0] # The center is always assumed to be 0 as of now
            V=rs[1]
            P=rs[2]

            sv=V.shape[0]
            aS=V.shape[1]

            Vp=np.matmul(self.computeUncertainMat()-self.Ac,V)

            U=np.zeros((self.n,1),dtype=object)

            for i in range(sv):
                s=0
                for j in range(aS):
                    #print(P[j][0],P[j][1])
                    s=s+(mp.mpi(P[j][0],P[j][1])*Vp[i][j])
                s=s+C[i]
                s_min=float(mp.nstr(s).split(',')[0][1:])
                s_max=float(mp.nstr(s).split(',')[1][:-1])
                U[i][0]=(s_min,s_max)

            #print(U)

            # Find a new Star set that overapproximates rs

            #print("Construction of the new overapproximated star is under construction!")



            semiDefFlag=False
            modelNewMax = Model("qp")
            modelNewMax.setParam('OutputFlag', False )

            C_new=np.zeros(self.n)
            V_new=V_randList[e]
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
                status = modelNewMax.Status
                if status==GRB.Status.UNBOUNDED:
                    0;
                    #print("UNBOUNDED")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        k=0
                        for v in modelNewMax.getVars():
                            #print('%s %g' % (v.varName, v.x))
                            UMax[k][0]=v.x
                            k=k+1

            except:
                print("Err")
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
                status = modelNewMin.Status
                if status==GRB.Status.UNBOUNDED:
                    0;
                    #print("UNBOUNDED")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        k=0
                        for v in modelNewMin.getVars():
                            #print('%s %g' % (v.varName, v.x))
                            UMin[k][0]=v.x
                            k=k+1

            except:
                print("Err")
                print("--")
                semiDefFlag=True


            ## DEBUG
            '''for d in range(UMin.shape[0]):
                if (UMin[d][0]>UMax[i][0]):
                    print(UMin)
                    print(UMax)
                    modelNewMin.write("minmodel.lp")
                    modelNewMax.write("maxmodel.lp")
                    exit(0)'''
            ##

            for i in range(self.n):
                #print((UMin[i][0],UMax[i][0]))
                #P_new.append((UMin[i][0],UMax[i][0]))
                P_new.append((min(UMin[i][0],UMax[i][0]),max(UMin[i][0],UMax[i][0])))

            #print(P_new)

            starNew=(C_new,V_new,P_new)

            '''print("-------Given Over-approximated Star-------")
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
            print("-----------------------")'''

            #exit(0)
            star_list.append(starNew)

        return star_list
        #-------------------------------------

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

    def computeUI_IntervalOld(self,rs):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2]


        sv=V.shape[0]
        aS=V.shape[1]

        Vp=np.matmul(self.computeUncertainMat()-self.Ac,V)

        U=np.zeros((self.n,1),dtype=object)

        for i in range(sv):
            s=0
            for j in range(aS):
                s=s+(mp.mpi(P[j][0],P[j][1])*Vp[i][j])
            s=s+C[i]
            s_min=float(mp.nstr(s).split(',')[0][1:])
            s_max=float(mp.nstr(s).split(',')[1][:-1])
            U[i][0]=(s_min,s_max)

        #print(U)

        # Find a new Star set that overapproximates rs
        semiDefFlag=False
        modelNewMax = Model("qp")
        modelNewMax.setParam('OutputFlag', False )

        #print("Construction of the new overapproximated star is under construction!")
        C_new=np.zeros(self.n)
        V_new=np.identity(self.n)
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
            #modelNewMax.write('dumpU.lp')
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
            #modelNewMin.write('dump.lp')
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

        # Manual Assert
        for i in range(self.n):
            if U[i][0][0]!=P_new[i][0] or U[i][0][1]!=P_new[i][1]:
                print("Found!")
                print(P_new)
                print(U)
                modelNewMin.write('dump.lp')
                exit(0)
        #--------------

        starNew=(C_new,V_new,P_new)

        '''print("-------Given Over-approximated Star-------")
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
        print("-----------------------")'''

        #exit(0)
        return starNew
        #-------------------------------------

    def computeUI_Interval(self,rs):
        '''
        This method aprroximates a reachable set with uncertainties represented
        as a star to another bloated star
        '''

        start_time=time.time()
        C=rs[0] # The center is always assumed to be 0 as of now
        V=rs[1]
        P=rs[2]


        sv=V.shape[0]
        aS=V.shape[1]

        Vp=np.matmul(self.computeUncertainMat()-self.Ac,V)

        U=np.zeros((self.n,1),dtype=object)

        for i in range(sv):
            s=0
            for j in range(aS):
                s=s+(mp.mpi(P[j][0],P[j][1])*Vp[i][j])
            s=s+C[i]
            s_min=float(mp.nstr(s).split(',')[0][1:])
            s_max=float(mp.nstr(s).split(',')[1][:-1])
            U[i][0]=(s_min,s_max)

        #print(U)
        C_new=np.zeros(self.n)
        V_new=np.identity(self.n)
        P_new=[]

        for i in range(self.n):
            #print((UMin[i][0],UMax[i][0]))
            #P_new.append((UMin[i][0],UMax[i][0]))
            P_new.append((U[i][0][0],U[i][0][1]))


        Profiling.comp_U=Profiling.comp_U+(time.time()-start_time)
        starNew=(C_new,V_new,P_new)

        '''print("-------Given Over-approximated Star-------")
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
        print("-----------------------")'''

        #exit(0)
        return starNew
        #-------------------------------------


    @staticmethod
    def joinBasisVecs(v1,v2):
        r=v1.shape[0]
        c=v1.shape[1]
        c2=v2.shape[1]
        V=np.zeros((r,c+c2))
        for i in range(r):
            for j in range(c):
                V[i][j]=v1[i][j]
        for i in range(r):
            j2=0
            for j in range(c,c+c2):
                V[i][j]=v2[i][j2]
                j2=j2+1
        return V

    @staticmethod
    def addStars(st1,st2):
        '''
        Given two stars st1 and st2, this functions performs
        minkowski sum of the two stars.
        '''

        start_time=time.time()
        C1=st1[0]
        C2=st2[0]
        C_new=list(map(add, C1, C2))
        #C_new=C1+C2

        V1=st1[1]
        V2=st2[1]
        V_new=CompU.joinBasisVecs(V1,V2)

        P1=st1[2]
        P2=st2[2]
        P_new=P1+P2

        Profiling.add_stars=Profiling.add_stars+(time.time()-start_time)

        #print("V_shape: ",V1.shape[1],V2.shape[1])
        #print("Ps: ",len(P1),len(P2))
        #print("P_new: ",len(P_new))
        #P_new=Visualization.landStars(P1,P2)

        return (C_new,V_new,P_new)

    @staticmethod
    def addStarsList(st1,st2):
        '''
        Given two stars st1 and st2, this functions performs
        minkowski sum of the two stars.
        '''
        l=[]

        for i in range(len(st1)):
            l.append(CompU.addStars(st1[i],st2[i]))

        return l

    @staticmethod
    def addStarsPred(st1,st2):
        '''
        Given two stars st1 and st2, this functions performs
        minkowski sum of the two stars.
        '''

        C1=st1[0]
        C2=st2[0]
        C_new=list(map(add, C1, C2))
        #C_new=C1+C2

        V1=st1[1]
        V2=st2[1]
        V_new=CompU.joinBasisVecs(V1,V2)

        P1=st1[2]
        P2=st2[2]
        P_new=CompU.andPred(P1,P2)

        #print("V_shape: ",V1.shape[1],V2.shape[1])
        #print("Ps: ",len(P1),len(P2))
        #print("P_new: ",len(P_new))
        #P_new=Visualization.landStars(P1,P2)

        return (C_new,V_new,P_new)

    @staticmethod
    def andPred(P1,P2):
        '''
        Given two predicates P1 and P2, this functions performs
        conjunction of two stars.
        '''

        pred1=P1[0]
        d1=P1[1]
        n1=pred1.shape[0]
        s1=pred1.shape[1]


        pred2=P2[0]
        d2=P2[1]
        n2=pred2.shape[0]
        s2=pred2.shape[1]

        P_new=np.zeros((n1+n2,s1+s2))
        D_new=np.zeros((n1+n2,1),dtype=object)

        for i in range(n1+n2):
            for j in range(s1+s2):
                if (i<n1 and j<s1):
                    P_new[i][j]=pred1[i][j]
                elif (i>n1 and j>s1):
                    P_new[i][j]=pred2[n1-i][s1-j]
            if i<n1:
                D_new[i][0]=d1[i][0]
            else:
                D_new[i][0]=d2[n1-i][0]

        return (P_new, D_new)

    @staticmethod
    def prodMatStars(M,RS):
        '''
        Given a matrix M and a star RS, perform M times RS
        '''

        start_time=time.time()
        C=RS[0]
        C_new=np.matmul(M,C)
        #C_new=C1+C2

        V=RS[1]
        V_new=np.matmul(M,V)

        P=RS[2]
        P_new=P
        #P_new=Visualization.landStars(P1,P2)
        Profiling.prod_mat_stars=Profiling.prod_mat_stars+(time.time()-start_time)
        #Profiling.prod_mat_stars=Profiling.prod_mat_stars+1

        return (C_new,V_new,P_new)

    @staticmethod
    def prodMatStarsList(M,RS_list):
        '''
        Given a matrix M and a star RS, perform M times RS
        '''

        l=[]

        for rs in RS_list:
            l.append(CompU.prodMatStars(M,rs))

        return l


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

    C2=[0,0]
    V2=np.array([
    [1,0],
    [0,1]
    ])
    P2=np.array([
    [0,-1],
    [1,1]
    ])
    D2=np.array([
    [(-1,1)],
    [(-1,1)]
    ])
    rs2=(C2,V2,(P2,D2))

    A=CompU.addStarsPred(rs,rs2)
    print(A[0])
    print(A[1])
    print(A[2][0])
    print(A[2][1])
    exit(0)

    pnew=np.array([
    [1,1],
    [1,-1]
    ])

    u=CompU(A,E)
    U=u.computeUI_Pred(rs,pnew)
    print(U[0])
    print(U[1])
    print(U[2][0])
    print(U[2][1])


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
    P=[(1,1),(1,1)]
    rs=(C,V,P)


    u=CompU(A,E)
    u.computeUI_Interval(rs)
    exit(0)

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
