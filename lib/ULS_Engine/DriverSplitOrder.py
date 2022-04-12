'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

import sys

from Benchmarks import *
from Consolidated import *
from SplitMet import *

mode='.'
h=0.001

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
        T=1500
        th1=1
        th2=2
        start=0
        n=200
        step=0.01
        methodList=['Kagstrom1','Kagstrom2','Loan']
        name="Stable1"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=(ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=1500
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="Stable2"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        #sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=1500
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="Stable3"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=2000
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="Stable4"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=50
        n=1000
        step=0.01
        name="FlightEnvelope"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")
        sp.printReachableSetOrdComp(th1,th2,ErB,name)

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
        T=200
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="Co-OpI"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")
        sp.printReachableSetOrdComp(th1,th2,ErB,name)

    def pkpd2():
        dynA=Benchmarks.PKPD2.A
        dynB=Benchmarks.PKPD2.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=200
        E={
        (0,4): [1-(P/100),1+(P/100)],
        #(2,0): [1-(P/100),1+(P/100)],
        #(2,2): [1-(P/100),1+(P/100)]
        (3,0): [1-(P/100),1+(P/100)],
        (3,3): [1-(P/100),1+(P/100)]
        }
        C=[5,5,5,5,5]
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
        T=200
        th1=1
        th2=2
        start=0
        n=50
        step=0.01
        name="PKPD2"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        #sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=2
        start=10
        n=20
        step=0.01
        name="DcConv"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        #o=OrdUnc(A)
        #(ErT,ErB)=o.printReportCompare()
        #sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        #sp.printReachableSetOrdComp(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="Space"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        #sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

    def holesCXc():
        dynA=Benchmarks.HolesCXc.A
        dynB=Benchmarks.HolesCXc.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=50
        E={
        (0,3): [1-(P/100),1+(P/100)],
        (1,2): [1-(P/100),1+(P/100)],
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
        T=100
        th1=0
        th2=1
        start=50
        n=2000
        step=0.01
        name="Holes"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        #o=OrdUnc(A)
        #(ErT,ErB)=o.printReportCompare()
        #sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        #sp.printReachableSetOrdComp(th1,th2,ErB,name)
        flow=([-1.115079e+00,1.120251e+00,1.120251e+00,-1.115079e+00,-1.115079e+00] , [-9.528530e-01,-9.528530e-01,1.142746e+00,1.142746e+00,-9.528530e-01])
        sp2=Split(A,E,rs,T)
        sp2.printReachableSetAll(th1,th2,flow,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=50
        n=5000
        step=0.01
        name="Motor1"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=4
        start=0
        n=200
        step=0.01
        name="Motor2"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="AircraftDynamics"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")
        sp.printReachableSetOrdComp(th1,th2,ErB,name)

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
        T=2000
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="GiradI-Full"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")
        sp.printReachableSetOrdComp(th1,th2,ErB,name)

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
        T=200
        th1=0
        th2=3
        start=0
        n=200
        step=0.01
        name="GiradII"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")
        sp.printReachableSetOrdComp(th1,th2,ErB,name)

    def acc():
        dynA=Benchmarks.ACC.A
        dynB=Benchmarks.ACC.B
        A=DriverCompU.createMatrix(dynA,dynB,mode,h)
        P=25
        E=Benchmarks.ACC.E
        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1],
        ])
        #P=[(0,35),(5,50),(0,35),(1,1)]
        P=[(1,11),(1,11),(1,11),(1,1)]
        rs=(C,V,P)
        T=100
        th1=1
        th2=2
        start=50
        n=2000
        step=0.01
        name="ACC"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        #o=OrdUnc(A)
        #(ErT,ErB)=o.printReportCompare()
        #sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        flow=([7.010404e-01,1.105693e+01,1.105693e+01,7.010404e-01,7.010404e-01] , [-1.284642e-02,-1.284642e-02,1.202995e+01,1.202995e+01,-1.284642e-02])
        sp2.printReachableSetAll(th1,th2,flow,name+"_Org")
        #sp.printReachableSetOrdComp(th1,th2,ErB,name)

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
        T=100
        th1=0
        th2=1
        start=0
        n=3
        step=0.01
        name="LaneChange"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        #o=OrdUnc(A)
        #(ErT,ErB)=o.printReportCompare()
        #sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        flow=([1.964808e+00,5.301050e+01,5.301050e+01,1.964808e+00,1.964808e+00] , [2.824269e+00,2.824269e+00,3.675731e+00,3.675731e+00,2.824269e+00])
        sp2.printReachableSetAll(th1,th2,flow,name+"_Org")
        #sp.printReachableSetOrdComp(th1,th2,ErB,name)

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
        T=200
        th1=0
        th2=1
        start=0
        n=200
        step=0.01
        name="CoOpII"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=0
        n=50
        step=0.01
        name="Five Vehicle"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")

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
        T=200
        th1=0
        th2=1
        start=0
        n=50
        step=0.01
        name="Mathias"
        print(">>>>>>>>>",name,"<<<<<<<<<\n\n")
        o=OrdUnc(A)
        (ErT,ErB)=o.printReportCompare()
        sp=Split(A,ErT,rs,T)
        #sp.printReachableSetTwo(th1,th2,ErB,name)
        sp.printReachableSetOrdComp(th1,th2,ErB,name)
        sp2=Split(A,E,rs,T)
        #sp2.printReachableSetAll(th1,th2,name+"_Org")



#Batch 1
#DriverCompU.stableSystem1()
#DriverCompU.stableSystem2()
#DriverCompU.stableSystem3()
#DriverCompU.stableSystem4()
#DriverCompU.dcConv()
#DriverCompU.giradI()
#DriverCompU.spaceCraftRndzvs()
#DriverCompU.giradII()
#DriverCompU.aircraftDynamics()
#DriverCompU.motorTransmission2()



#Batch 2
#DriverCompU.acc()
#DriverCompU.flightEnvelope()
#DriverCompU.coOPVehiclesII()



#Batch 3
DriverCompU.laneChange()
#DriverCompU.fiveVehiclePlatton()



#Batch 4
#DriverCompU.holesCXc()
#DriverCompU.pkpd2()
#DriverCompU.motorTransmission1()
#DriverCompU.coOPVehiclesI()



#Batch for Zono (Ac Basis)
#DriverCompU.stableSystem1()
#DriverCompU.stableSystem2()
#DriverCompU.stableSystem3()
#DriverCompU.stableSystem4()
#DriverCompU.dcConv()
#DriverCompU.giradI()
#DriverCompU.spaceCraftRndzvs()
#DriverCompU.giradII()
#DriverCompU.aircraftDynamics()
#DriverCompU.motorTransmission2()
#DriverCompU.acc()
#DriverCompU.flightEnvelope()
#DriverCompU.coOPVehiclesII()
#DriverCompU.laneChange()
#DriverCompU.fiveVehiclePlatton()
#DriverCompU.holesCXc()
#DriverCompU.pkpd2()
#DriverCompU.motorTransmission1()
#DriverCompU.mathias()


#DriverCompU.holesCXc()
