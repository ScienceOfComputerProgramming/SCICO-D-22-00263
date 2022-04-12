'''
This file provides APIs to drive the Robust Metric
computation.
'''
from RobustMetric import *
from Benchmarks import *
import time

mode='.'
h=0.01

class DriverRM:
    '''
    APIs to compute the robustness metric for various Benchmarks
    '''
    @staticmethod
    def createMatrix(A,B,mode,h):
        ''' Creates a single matrix based on
        . or +.
        In case of . a rough approximation is
        done'''

        n1=np.size(A,0)
        if (np.size(B)>0):
            n2=np.size(B,1)
        else:
            n2=0
        n=n1+n2

        C=np.zeros((n,n),dtype=np.float)
        if mode=='+':
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1
        elif mode=='.':
            I=np.zeros((n1,n1),dtype=np.float)
            for i in range(n1):
                I[i][i]=1
            A2=h*A
            A2=np.add(I,A2)
            B2=h*B
            for i in range(n1):
                for j in range(n1):
                    C[i][j]=A2[i][j]
            for i in range(n1):
                j2=0
                for j in range(n1,n1+n2):
                    C[i][j]=B2[i][j2]
                    j2=j2+1
            for i in range(n1,n1+n2):
                C[i][i]=1

        return C

    def dcConv():
        dynA=Benchmarks.DCConv.A
        dynB=Benchmarks.DCConv.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,2),(1,1)]
        #cells=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1]
        ])
        P=[(-1,1),(-1,1),(1,1)]
        rs=(C,V,P)
        T=100
        th1=0
        th2=2
        name="DcConv"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(3,5),(-1,1),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def holesCXc():
        dynA=Benchmarks.HolesCXc.A
        dynB=Benchmarks.HolesCXc.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,3),(1,2),(3,2),(4,3),(5,5),(9,9)]
        '''cells=[]
        for i in range(10):
            for j in range(10):
                cells.append((i,j))'''
        C=[0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1],
        ])
        P=[(-1,1),(-1,1),(1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V,P)
        T=100
        th1=0
        th2=1
        name="HolesCXc"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.05,1.5),(-1,1),(1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def acc():
        dynA=Benchmarks.ACC.A
        dynB=Benchmarks.ACC.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        P=10
        cells=[(0,3),(2,3)]
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        P=[(0,35),(5,50),(0,35),(1,1)]
        rs=(C,V,P)
        T=100
        th1=0
        th2=1
        name="ACC"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(42,47),(5,50),(0,35),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def laneChange():
        dynA=Benchmarks.LaneChange.A
        dynB=Benchmarks.LaneChange.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(3,6),(4,6),(5,6)]
        C=[0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0],
        [0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1],
        ])
        P=[(0,50),(3,3.5),(0,0),(20,30),(0,0),(0,0),(1,1)]
        rs=(C,V,P)
        T=100
        th1=0
        th2=1
        start=0
        n=3
        step=0.01
        name="LaneChange"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(53,60),(3,3.5),(0,0),(20,30),(0,0),(0,0),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def pkpd2():
        dynA=Benchmarks.PKPD2.A
        dynB=Benchmarks.PKPD2.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,4),(0,0)]
        C=[5,5,5,5,5]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        rs=(C,V,P)
        T=100
        name="PKPD2"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(8,9),(0,10),(0,10),(1,8),(0,200)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def giradI():
        dynA=Benchmarks.GiradI.A
        dynB=Benchmarks.GiradI.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,0),(1,0)]
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(0.9,1.1),(-0.1,0.1)]
        rs=(C,V,P)
        T=10
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="GiradI-Full"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.2,1.7),(-0.1,0.1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def giradII():
        dynA=Benchmarks.GiradII.A
        dynB=Benchmarks.GiradII.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,0),(2,2),(4,4)]
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(0.9,1.1),(-0.1,0.1),(0.9,1.1),(-0.1,0.1),(0.9,1.1)]
        rs=(C,V,P)
        T=100
        name="GiradII"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.2,1.7),(-0.1,0.1),(0.9,1.1),(-0.1,0.1),(0.9,1.1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def fiveVehiclePlatton():
        dynA=Benchmarks.FiveVehiclePlatton.A
        dynB=Benchmarks.FiveVehiclePlatton.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (3,7): [1-(P/100),1+(P/100)],
        (4,6): [1-(P/100),1+(P/100)]
        }
        E2={
        (2,1): [1-(P/100),1+(P/100)],
        (5,1): [1-(P/100),1+(P/100)],
        (5,3): [1-(P/100),1+(P/100)],
        (8,3): [1-(P/100),1+(P/100)],
        (3,7): [1-(P/100),1+(P/100)],
        (4,6): [1-(P/100),1+(P/100)]
        }
        cells=[(2,1),(5,1),(5,3),(8,3),(3,7),(4,6)]
        C=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        VI=np.identity(15)
        P=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)]
        rs=(C,VI,P)
        T=100
        name="5-Veh"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.01,1.2),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)]
        unsafe=(C,VI,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def motorTransmission1():
        dynA=Benchmarks.MotorTransmission1.A
        dynB=Benchmarks.MotorTransmission1.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,6),(1,6)]
        C=[0,0,0,0,0,0,0]
        V2=np.array([
        [1,0,0,0,0,0,0],
        [0,1,0,0,0,0,0],
        [0,0,1,0,0,0,0],
        [0,0,0,1,0,0,0],
        [0,0,0,0,1,0,0],
        [0,0,0,0,0,1,0],
        [0,0,0,0,0,0,1]
        ])
        P=[(-1,1),(-1,1),(-1,-1),(-1,1),(-1,1),(1,1),(1,1)]
        rs=(C,V2,P)
        T=100
        name="Motor 1"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=P=[(1.005,2),(-1,1),(-1,-1),(-1,1),(-1,1),(1,1),(1,1)]
        unsafe=(C,V2,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def spaceCraftRndzvs():
        dynA=Benchmarks.SpaceCraftRndzvs.A
        dynB=Benchmarks.SpaceCraftRndzvs.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(2,1),(2,0),(2,3),(3,2)]
        C=[0,0,0,0,0,0]
        V2=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1],
        ])
        P=[(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V2,P)
        T=100
        name="Space"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.1,1.5),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        unsafe=(C,V2,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def aircraftDynamics():
        dynA=Benchmarks.AircraftDynamics.A
        dynB=Benchmarks.AircraftDynamics.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(2,3),(3,2)]
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(-1,1),(-1,1),(20,30),(20,30)]
        rs=(C,V,P)
        T=100
        name="AircraftDynamics"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(-1,1),(-1,1),(32,34),(20,30)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def coOPVehiclesI():
        dynA=Benchmarks.CoOPVehiclesI.A
        dynB=Benchmarks.CoOPVehiclesI.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(0,0),(2,0),(2,1),(2,8),(5,0),(5,8)]
        C=[0,0,0,0,0,0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,1]
        ])
        P=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(1,1)]
        rs=(C,V,P)
        T=100
        name="Co-Op I"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(1.3,1.5),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(1,1)]
        unsafe=(C,V,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)

    def flightEnvelope():
        dynA=Benchmarks.FlightEnvelope.A
        dynB=Benchmarks.FlightEnvelope.B
        A=DriverRM.createMatrix(dynA,dynB,mode,h)
        cells=[(3,7),(4,6)]
        C=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        VI=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        ])
        T=100
        P=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,VI,P)
        name="FlightEnvelope"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        P_unsafe=[(-1,1),(-1,1),(-1,1),(1.2,1.4),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        unsafe=(C,VI,P_unsafe)
        rm=RobustMetric(A,cells)
        time_taken=time.time()
        metric=rm.getRobustMetricAll(rs,T,unsafe)
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        print(metric)
        print("Time Taken: ",time.time()-time_taken)




DriverRM.motorTransmission1()
