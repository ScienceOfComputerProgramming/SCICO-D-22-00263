'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

import sys

from Benchmarks import *
from Consolidated import *

mode='.'
h=0.01

class DriverCompU:

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

    def stableSystem1():
        dynA=Benchmarks.StableSystem1.A
        dynB=Benchmarks.StableSystem1.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (2,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0]
        V=np.array([
        [1,-1,1],
        [-1,1,0],
        [1,0,-1],
        ])
        P=[(-1,1),(-1,1),(-1,1)]
        rs=(C,V,P)
        T=2050
        th1=1
        th2=2
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Stable1")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def stableSystem2():
        dynA=Benchmarks.StableSystem2.A
        dynB=Benchmarks.StableSystem2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(-1,1),(-1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Stable2")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def stableSystem3():
        dynA=Benchmarks.StableSystem3.A
        dynB=Benchmarks.StableSystem3.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(-1,1),(-1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Stable3")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def stableSystem4():
        dynA=Benchmarks.StableSystem4.A
        dynB=Benchmarks.StableSystem4.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0]
        V=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1],
        ])
        P=[(-1,1),(-1,1),(1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Stable4")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def flightEnvelope():
        dynA=Benchmarks.FlightEnvelope.A
        dynB=Benchmarks.FlightEnvelope.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (3,7): [1-(P/100),1+(P/100)],
        (4,6): [1-(P/100),1+(P/100)]
        }
        E2={
        (0,1): [1-(P/100),1+(P/100)],
        (0,3): [1-(P/100),1+(P/100)],
        (1,4): [1-(P/100),1+(P/100)],
        (8,11): [1-(P/100),1+(P/100)]
        }
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
        V=np.array([
        [1,0,0,0,1,0,0,1,0,0,0,-1,0,0,1,1],
        [0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0],
        [-1,0,1,0,0,0,-1,0,1,0,1,0,0,-1,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,-1,0,0,1,0,0,0,1,0,-1,0,1,0,0,0],
        [1,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,-1,0,0,0,0,-1,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,-1,0,0,0,0,0,0,1,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0],
        [0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,-1],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,-1,0,0,0,0,-1,0,0,0,0,0,0,1,0,0],
        [0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,1,0],
        [1,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,1]
        ])
        P=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,VI,P)
        T=2050
        th1=0
        th2=1
        start=50
        n=1000
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"FlightEnvelope20p")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def coOPVehiclesI():
        dynA=Benchmarks.CoOPVehiclesI.A
        dynB=Benchmarks.CoOPVehiclesI.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,0): [1-(P/100),1+(P/100)],
        (2,1): [1-(P/100),1+(P/100)],
        (2,8): [1-(P/100),1+(P/100)],
        (5,0): [1-(P/100),1+(P/100)],
        (5,8): [1-(P/100),1+(P/100)]
        }
        E2={
        (2,0): [1-(P/100),1+(P/100)],
        (5,0): [1-(P/100),1+(P/100)]
        }
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
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Loan']
        co=SplitBloat(E,rs,T,"CoOpI")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def pkpd2():
        dynA=Benchmarks.PKPD2.A
        dynB=Benchmarks.PKPD2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (0,4): [1-(P/100),1+(P/100)],
        #(2,0): [1-(P/100),1+(P/100)],
        #(2,2): [1-(P/100),1+(P/100)]
        (3,0): [1-(P/100),1+(P/100)],
        (3,3): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        V2=np.array([
        [1,0,-1,0,-1],
        [0,1,1,0,0],
        [1,1,0,0,-1],
        [1,0,-1,1,0],
        [1,1,0,1,1],
        ])
        P=[(1,6),(0,10),(0,10),(1,8),(0,200)]
        rs=(C,V,P)
        T=2050
        th1=1
        th2=2
        start=50
        n=500
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"PKPD2-3")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def dcConv():
        dynA=Benchmarks.DCConv.A
        dynB=Benchmarks.DCConv.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E={
        (0,2): [1-(P/100),1+(P/100)],
        (1,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0]
        V2=np.array([
        [1,0,0],
        [0,1,0],
        [0,0,1]
        ])
        V=np.array([
        [1,-1,0],
        [0,1,-1],
        [-1,0,1]
        ])
        P=[(-1,1),(-1,1),(1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=2
        start=10
        n=20
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"DC_10p")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def spaceCraftRndzvs():
        dynA=Benchmarks.SpaceCraftRndzvs.A
        dynB=Benchmarks.SpaceCraftRndzvs.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (2,1): [1-(P/100),1+(P/100)],
        (2,0): [1-(P/100),1+(P/100)],
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)],
        }
        C=[0,0,0,0,0,0]
        V2=np.array([
        [1,0,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,1,0,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1],
        ])
        V=np.array([
        [1,0,0,1,0,-1],
        [0,-1,0,1,0,-1],
        [-1,-1,1,0,0,1],
        [0,1,0,1,1,-1],
        [-1,-1,0,0,1,0],
        [1,-1,0,-1,0,1],
        ])
        P=[(-1,1),(-1,1),(1,1),(1,1),(1,1),(1,1)]
        rs=(C,V2,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Spacev2")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def holesCXc():
        dynA=Benchmarks.HolesCXc.A
        dynB=Benchmarks.HolesCXc.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (0,3): [1-(P/100),1+(P/100)],
        (1,2): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)],
        (4,3): [1-(P/100),1+(P/100)]
        }
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
        T=2050
        th1=0
        th2=1
        start=50
        n=2000
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Holes20p")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def motorTransmission1():
        dynA=Benchmarks.MotorTransmission1.A
        dynB=Benchmarks.MotorTransmission1.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (0,6): [1-(P/100),1+(P/100)],
        (1,6): [1-(P/100),1+(P/100)]
        }
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
        V=np.array([
        [1,0,0,1,0,0,-1],
        [0,1,0,-1,0,-1,0],
        [1,0,1,1,1,-1,0],
        [0,-1,0,1,0,0,-1],
        [-1,0,0,0,1,0,1],
        [0,0,-1,0,0,1,0],
        [1,0,0,1,0,0,1]
        ])
        P=[(-1,1),(-1,1),(-1,-1),(-1,1),(-1,1),(1,1),(1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=50
        n=5000
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Motor1_20p")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def motorTransmission2():
        dynA=Benchmarks.MotorTransmission2.A
        dynB=Benchmarks.MotorTransmission2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (4,0): [1-(P/100),1+(P/100)],
        (4,1): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1]
        ])
        P=[(-1,1),(-1,1),(1,1),(1,1),(-1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=4
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Motor2")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def aircraftDynamics():
        dynA=Benchmarks.AircraftDynamics.A
        dynB=Benchmarks.AircraftDynamics.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=80
        E={
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)]
        }
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(-1,1),(-1,1),(20,30),(20,30)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"AircraftDynamics")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def giradI():
        dynA=Benchmarks.GiradI.A
        dynB=Benchmarks.GiradI.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (1,0): [1-(P/100),1+(P/100)]
        }
        C=[0,0]
        V=np.array([
        [1,0],
        [0,1]
        ])
        P=[(0.9,1.1),(-0.1,0.1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"GiradI")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def giradII():
        dynA=Benchmarks.GiradII.A
        dynB=Benchmarks.GiradII.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=2
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,2): [1-(P/100),1+(P/100)],
        (4,4): [1-(P/100),1+(P/100)]
        }
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
        T=2050
        th1=0
        th2=3
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"GiradII")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def acc():
        dynA=Benchmarks.ACC.A
        dynB=Benchmarks.ACC.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E=Benchmarks.ACC.E
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        P=[(0,35),(5,50),(0,35),(1,1)]
        rs=(C,V,P)
        T=2050
        th1=1
        th2=2
        start=50
        n=2000
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"ACC")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'slow')

    def laneChange():
        dynA=Benchmarks.LaneChange.A
        dynB=Benchmarks.LaneChange.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=10
        E=Benchmarks.LaneChange.E
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
        T=2050
        th1=0
        th2=1
        start=50
        n=500
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"LaneChange")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def coOPVehiclesII():
        dynA=Benchmarks.CoOPVehiclesII.A
        dynB=Benchmarks.CoOPVehiclesII.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E={
        (0,0): [1-(P/100),1+(P/100)],
        (2,0): [1-(P/100),1+(P/100)],
        (2,1): [1-(P/100),1+(P/100)],
        (2,8): [1-(P/100),1+(P/100)],
        (5,0): [1-(P/100),1+(P/100)],
        (5,8): [1-(P/100),1+(P/100)]
        }
        E2={
        (2,0): [1-(P/100),1+(P/100)],
        (5,0): [1-(P/100),1+(P/100)]
        }
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
        P=[(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1),(-1,1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"CoOpII")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def fiveVehiclePlatton():
        dynA=Benchmarks.FiveVehiclePlatton.A
        dynB=Benchmarks.FiveVehiclePlatton.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
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
        C=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        VI=np.array([
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
        ])
        V=np.array([
        [1,0,0,0,1,0,0,1,0,0,0,-1,0,0,1,1],
        [0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0],
        [-1,0,1,0,0,0,-1,0,1,0,1,0,0,-1,0,0],
        [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,-1,0,0,1,0,0,0,1,0,-1,0,1,0,0,0],
        [1,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,1,0,0,-1,0,0,0,0,-1,0],
        [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
        [0,-1,0,0,0,0,0,0,1,0,0,0,0,1,0,0],
        [0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0],
        [0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
        [1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,-1],
        [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
        [0,-1,0,0,0,0,-1,0,0,0,0,0,0,1,0,0],
        [0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,1,0],
        [1,0,0,0,0,0,-1,0,0,0,0,0,-1,0,0,1]
        ])
        P=[(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1),(0.9,1.1)]
        rs=(C,VI,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=50
        step=0.01
        methodList=['Kagstrom1','Loan']
        #sp=Split(A,E2,rs,T)
        #sp.printReachableSetAll(0,1,"FiveVehiclePlatton")
        co=SplitBloat(E2,rs,T,"FiveVehiclePlatton")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')

    def mathias():
        dynA=Benchmarks.Mathias.A
        dynB=Benchmarks.Mathias.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=20
        E=Benchmarks.Mathias.E
        C=[0,0,0,0,0]
        V=np.array([
        [1,0,0,0,0],
        [0,1,0,0,0],
        [0,0,1,0,0],
        [0,0,0,1,0],
        [0,0,0,0,1],
        ])
        P=[(-0.1,0.1),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1),(-0.1,0.1)]
        rs=(C,V,P)
        T=2050
        th1=0
        th2=1
        start=0
        n=50
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        co=SplitBloat(E,rs,T,"Mathias")
        co.getStats(dynA,dynB,th1,th2,start,n,step,methodList,'fast')


'''
#Batch 1
DriverCompU.stableSystem1()
DriverCompU.stableSystem2()
DriverCompU.stableSystem3()
DriverCompU.stableSystem4()
#DriverCompU.dcConv()
DriverCompU.giradI()
DriverCompU.spaceCraftRndzvs()
DriverCompU.giradII()
DriverCompU.aircraftDynamics()
DriverCompU.motorTransmission2()
'''

'''
#Batch 2
DriverCompU.acc()
DriverCompU.flightEnvelope()
DriverCompU.coOPVehiclesII()
'''

'''
#Batch 3
DriverCompU.laneChange()
DriverCompU.fiveVehiclePlatton()
'''

'''
#Batch 4
DriverCompU.holesCXc()
DriverCompU.pkpd2()
DriverCompU.motorTransmission1()
DriverCompU.coOPVehiclesI()
'''


#Batch for Zono (Ac Basis)
DriverCompU.stableSystem1()
DriverCompU.stableSystem2()
DriverCompU.stableSystem3()
DriverCompU.stableSystem4()
DriverCompU.dcConv()
DriverCompU.giradI()
DriverCompU.spaceCraftRndzvs()
DriverCompU.giradII()
DriverCompU.aircraftDynamics()
DriverCompU.motorTransmission2()
DriverCompU.acc()
DriverCompU.flightEnvelope()
DriverCompU.coOPVehiclesI()
DriverCompU.laneChange()
DriverCompU.fiveVehiclePlatton()
DriverCompU.holesCXc()
DriverCompU.pkpd2()
DriverCompU.motorTransmission1()
DriverCompU.laneChange()


#DriverCompU.motorTransmission1()
