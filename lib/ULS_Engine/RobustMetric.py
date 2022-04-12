'''
This file provides APIs to compute the robustness metric of
a linear uncertain system, given an unsafe set
'''

from OrderUncertainties import *
from SplitMet import *
from VisualizationReachSet import *

from BloatAPI import *
from Consolidated import *

class RobustMetric:
    '''
    Compute the robustness metric of a given linear uncertain
    system and a corresponding unsafe set
    '''

    def __init__(self,A,cells):
        self.A=A # Linear system (without perturbation)
        self.cells=cells

    def distributePerturbation(self,p,choice=3):
        '''
        Given a perturbation budget `p`, disribute it in the cells
        mentioned in `self.cells`
        cells=[...,(i,j),...], where (i,j) is a cell in `self.A`

        choices:
            - 1: Linear Distribution
            - 2: Harmonic Distribution
            - 3: Equal Distribution
            - 4: Return All
        '''

        if choice==1:
            Er=self.distributePerturbationLinear(p)
            return Er
        elif choice==2:
            Er=self.distributePerturbationHarmonic(p)
            return Er
        elif choice==3:
            Er=self.distributePerturbationEqual(p)
            return Er
        elif choice==4:
            ErLin=self.distributePerturbationLinear(p)
            ErHar=self.distributePerturbationHarmonic(p)
            ErEq=self.distributePerturbationEqual(p)
            return (ErLin,ErHar,ErEq)


    def distributePerturbationLinear(self,p):
        '''
        Given a perturbation budget `p`, disribute it in the cells
        mentioned in `self.cells`
        cells=[...,(i,j),...], where (i,j) is a cell in `self.A`

        # Algorithm

        - for a cell c_i \in cells, plug p, B[c_i]=p; and call the ordering algorithm
        - obtain the \sigma_{c_i} from the ordering algoritm
        - order the cells based on the obtained sigma_ values from the above step
            - highest sigma, gets the lowest propotion of p, and so.
        '''
        Er={}
        sensitivityMat=OrdUnc(self.A).getSVSentivity(p)
        ep=np.amin(sensitivityMat)
        sm=0
        for i in self.cells:
            if sensitivityMat[i]==0:
                sensitivityMat[i]=1e-20
            sm+=sensitivityMat[i]
        sortedCells=self.sortCells(sensitivityMat)
        for i in range(len(sortedCells)):
            tmp=(p/sm)*sensitivityMat[sortedCells[-1-i]]
            #print(self.cells[i],": ",tmp)
            Er[sortedCells[i]]=[(100-tmp)/100,(100+tmp)/100]
        #print(Er)

        return Er

    def sortCells(self,sensitivityMat):
        '''
        sort cells based on sensitivityMat values
        '''

        n_cells=len(self.cells)
        cellsOrg=self.cells.copy()
        for i in range(n_cells):
            indx=i
            mx=-9999
            for j in range(i,n_cells):
                if sensitivityMat[cellsOrg[j]]>mx:
                    mx=sensitivityMat[cellsOrg[j]]
                    indx=j
            #print(cellsOrg[indx])
            tm=cellsOrg[i]
            cellsOrg[i]=cellsOrg[indx]
            cellsOrg[indx]=tm
        return cellsOrg

    def distributePerturbationHarmonic(self,p):
        '''
        Given a perturbation budget `p`, disribute it in the cells
        mentioned in `self.cells`
        cells=[...,(i,j),...], where (i,j) is a cell in `self.A`

        # Algorithm

        - for a cell c_i \in cells, plug p, B[c_i]=p; and call the ordering algorithm
        - obtain the \sigma_{c_i} from the ordering algorithm
        - each cell gets perturbation based on the proportion of 1/sigma_{c_i}
        '''
        Er={}
        sensitivityMat=OrdUnc(self.A).getSVSentivity(p)
        sm=0
        for i in self.cells:
            if sensitivityMat[i]==0:
                sensitivityMat[i]=1e-20
            sm+=float(1/sensitivityMat[i])

        for i in range(len(self.cells)):
            tmp=float(p/sm)*float(1/sensitivityMat[self.cells[i]])
            #print(self.cells[i],": ",tmp)
            Er[self.cells[i]]=[(100-tmp)/100,(100+tmp)/100]

        return Er

    def distributePerturbationEqual(self,p):
        '''
        Given a perturbation budget `p`, disribute it in the cells
        mentioned in `self.cells`
        cells=[...,(i,j),...], where (i,j) is a cell in `self.A`

        # Algorithm

        - all cells get equal share of perturbation
        '''
        Er={}
        n_cells=len(self.cells)
        for i in range(n_cells):
            tmp=p/n_cells
            #print(self.cells[i],": ",tmp)
            Er[self.cells[i]]=[(100-tmp)/100,(100+tmp)/100]
        #print(Er)

        return Er

    def checkIntersection(star,unsafe):
        '''
        Returns:
            - False, if `star` and `unsafe` doesn't intersect
            - True, otherwise

        Algorithm:
            - Encode `star` and `unsafe` as Gurobi variables
            - Find a point which lies on both `star` and `unsafe`
                - If such a point exists, i.e. the optimization is feasible: return True
                - Otherwise: return False
        '''
        intersectFlag=True

        model = Model("qp")

        # Encode `star` . . .

        C_star=star[0]
        V_star=star[1]
        P_star=star[2]
        n_star=len(C_star)
        vecSize_star=V_star.shape[1]

        ## Encode predicate variables
        predVars_star=[]
        for i in range(vecSize_star):
            name="pred_"+str(i)
            predVars_star.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
            model.addConstr(predVars_star[i]>=P_star[i][0],name+".1")
            model.addConstr(predVars_star[i]<=P_star[i][1],name+".2")

        ## Encode the states
        stateVars_star=[]
        for i in range(n_star):
            #stateVars_star.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name="State_RS_"+str(i),vtype='C'))
            obj=C_star[i]
            for j in range(vecSize_star):
                obj=obj+(V_star[i][j]*predVars_star[j])
            #model.addConstr(stateVars_star[i]==obj)
            stateVars_star.append(obj)

        # . . . Encode `star` end

        # Encode `unsafe` . . .

        C_unsafe=unsafe[0]
        V_unsafe=unsafe[1]
        P_unsafe=unsafe[2]
        n_unsafe=len(C_unsafe)
        vecSize_unsafe=V_unsafe.shape[1]

        ## Encode predicate variables
        predVars_unsafe=[]
        for i in range(vecSize_unsafe):
            name="predU_"+str(i)
            predVars_unsafe.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
            model.addConstr(predVars_unsafe[i]>=P_unsafe[i][0],name+".1")
            model.addConstr(predVars_unsafe[i]<=P_unsafe[i][1],name+".2")

        ## Encode the states
        stateVars_unsafe=[]
        for i in range(n_unsafe):
            #stateVars_unsafe.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name="State_Unsafe_"+str(i),vtype='C'))
            obj=C_unsafe[i]
            for j in range(vecSize_unsafe):
                obj=obj+(V_unsafe[i][j]*predVars_unsafe[j])
            #model.addConstr(stateVars_unsafe[i]==obj)
            stateVars_unsafe.append(obj)

        # . . . Encode `unsafe` end

        # Add constraints for finding a point that lies on both `star` and `unsafe`
        for i in range(n_star):
            model.addConstr(stateVars_star[i]==stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]>=stateVars_unsafe[i])
            #model.addConstr(stateVars_star[i]<=stateVars_unsafe[i])

        # Set Objective
        model.setObjective(stateVars_star[0]) # Can be anything --- we just care about feasibility

        # Check feasibility
        model.optimize()

        print()

        status = model.Status
        if status==GRB.Status.UNBOUNDED:
            print("UNBOUNDED ")
        else:
            if status == GRB.Status.INF_OR_UNBD or \
               status == GRB.Status.INFEASIBLE  or \
               status == GRB.Status.UNBOUNDED:
                intersectFlag=False
            else:
                intersectFlag=True
                '''for v in model.getVars():
                    print('%s %g' % (v.varName, v.x))'''

        return intersectFlag

    def getRobustMetric(self,IS,t,unsafe):
        '''
        Given initial set, unsafe set: return the robustness metric

        # Algorithm:
        - Keep increasing p:
            - \Lambda=self.distributePerturbation(p)
            - if unsafe:
                break
        '''
        p=0.1
        delta=2

        it=1

        t_taken=time.time()
        while True:
            #print(p)
            Er=self.distributePerturbation(p)
            #Er={}
            #print(Er)
            sp=Split(self.A,Er,IS,t)
            RS=sp.getReachableSetAll()
            fg_intersect=RobustMetric.checkIntersection(RS,unsafe)
            if fg_intersect==True:
                t_taken=time.time()-t_taken
                if True:
                    print("* Number of iterations: ",it)
                    print("* Total Perturbation: ",p/delta)
                    print("* Total Time Taken: ",t_taken)
                Visualization.displayRSnUnsafe(RS,unsafe,th1=0,th2=1,name="RS")
                return ((p/delta),Er)
            #Visualization.displayRSnUnsafe(RS,unsafe,th1=0,th2=1,name="RS")
            p*=delta
            it+=1

        return -404

    def compareHeus(self,p,IS,t):
        '''
        Compare the three heuristics
        '''
        (ErLin,ErHar,ErEq)=self.distributePerturbation(p,choice=4)
        sp=Split(self.A,ErLin,IS,t)
        RS_lin=sp.getReachableSetAll()
        sp2=Split(self.A,ErHar,IS,t)
        RS_har=sp2.getReachableSetAll()
        sp3=Split(self.A,ErEq,IS,t)
        RS_eq=sp3.getReachableSetAll()

        Visualization.displayCompHeu(RS_lin,RS_har,RS_eq,th1=0,th2=1,name="Heu")

    def getRobustMetricAll(self,IS,t,unsafe):
        '''
        Given initial set, unsafe set: return the robustness metric

        # Algorithm:
        - Keep increasing p:
            - \Lambda=self.distributePerturbation(p)
            - if unsafe:
                break
        '''
        p=10
        delta=3

        it=1

        t_taken=time.time()
        er_dict={'lin':[],'har':[],'eq':[]}
        fg_lin=True
        fg_har=True
        fg_eq=True
        while True:
            Er_lin=self.distributePerturbation(p,choice=1)
            Er_har=self.distributePerturbation(p,choice=2)
            Er_eq=self.distributePerturbation(p,choice=3)
            sp_lin=Split(self.A,Er_lin,IS,t)
            sp_har=Split(self.A,Er_har,IS,t)
            sp_eq=Split(self.A,Er_eq,IS,t)
            RS_lin=sp_lin.getReachableSetAll()
            RS_har=sp_har.getReachableSetAll()
            RS_eq=sp_eq.getReachableSetAll()
            fg_intersect_lin=RobustMetric.checkIntersection(RS_lin,unsafe)
            fg_intersect_har=RobustMetric.checkIntersection(RS_har,unsafe)
            fg_intersect_eq=RobustMetric.checkIntersection(RS_eq,unsafe)
            if fg_lin==True and fg_intersect_lin==True:
                matE=SplitBloat(Er_lin,[],0,"c").matrixify(self.A)
                nrm=BloatKagstrom(self.A,matE).intervalNorm(p='fast')
                er_dict['lin']=[(p/delta),nrm,Er_lin]
                fg_lin=False
            if fg_har==True and fg_intersect_har==True:
                matE=SplitBloat(Er_har,[],0,"c").matrixify(self.A)
                nrm=BloatKagstrom(self.A,matE).intervalNorm(p='fast')
                er_dict['har']=[(p/delta),nrm,Er_har]
                fg_har=False
            if fg_eq==True and fg_intersect_eq==True:
                matE=SplitBloat(Er_eq,[],0,"c").matrixify(self.A)
                nrm=BloatKagstrom(self.A,matE).intervalNorm(p='fast')
                er_dict['eq']=[(p/delta),nrm,Er_eq]
                fg_eq=False
            if fg_eq==False and fg_har==False and fg_lin==False:
                return er_dict
            p*=delta
            it+=1

        return -404




# Test Case
if False:
    A=np.array([
    [20,3],
    [2,51.2]
    ])
    cells=[(0,1),(1,0)]
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=[(-1,1),(-1,1)]
    rs=(C,V,P)
    t=2
    P_unsafe=[(600,800),(600,800)]
    unsafe=(C,V,P_unsafe)
    rm=RobustMetric(A,cells)
    #rm.distributePerturbation(10)
    metric=rm.getRobustMetric(rs,t,unsafe)
    print(metric)
