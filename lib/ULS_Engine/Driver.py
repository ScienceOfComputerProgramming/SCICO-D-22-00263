'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

from ReachSetAPI import *
from Benchmarks import *

class Driver:

    def stableSystem1():
        A=Benchmarks.StableSystem1.A
        B=Benchmarks.StableSystem1.B
        mode='.'
        E={
        (0,2): [0.9,1.1],
        (2,1): [0.9,1.1],
        }
        IS=np.array([
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Stable System 1 <<<<<<<<<<<<<\n")
        rs=Split
        rs.compareReachSets()
        print("\n==========================================\n")

    def stableSystem2():
        A=Benchmarks.StableSystem2.A
        B=Benchmarks.StableSystem2.B
        mode='.'
        E={
        (0,1): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Stable System 2 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def stableSystem3():
        A=Benchmarks.StableSystem3.A
        B=Benchmarks.StableSystem3.B
        mode='.'
        E={
        (0,1): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Stable System 3 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def stableSystem4():
        A=Benchmarks.StableSystem4.A
        B=Benchmarks.StableSystem4.B
        mode='.'
        E={
        (0,2): [0.9,1.1],
        (2,1): [0.9,1.1],
        }
        IS=np.array([
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Stable System 4 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def flightEnvelope():
        A=Benchmarks.FlightEnvelope.A
        B=Benchmarks.FlightEnvelope.B
        mode='.'
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        E={
        (3,7): [0.9,1.1],
        (4,6): [0.9,1.1]
        }
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Flight Envelope <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def coOPVehiclesI():
        A=Benchmarks.CoOPVehiclesI.A
        B=Benchmarks.CoOPVehiclesI.B
        mode='.'
        E={
        (2,0): [0.9,1.1],
        (2,1): [0.9,1.1],
        (2,8): [0.9,1.1],
        (5,0): [0.9,1.1],
        (5,8): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Co Op Vehicles I <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def coOPVehiclesII():
        A=Benchmarks.CoOPVehiclesII.A
        B=Benchmarks.CoOPVehiclesII.B
        mode='.'
        E={
        (2,0): [0.9,1.1],
        (2,1): [0.9,1.1],
        (2,8): [0.9,1.1],
        (5,0): [0.9,1.1],
        (5,8): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Co Op Vehicles II <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def pkpd2():
        A=Benchmarks.PKPD2.A
        B=Benchmarks.PKPD2.B
        mode='.'
        E={
        (0,4): [0.9,1.1],
        (3,3): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> PKPD 2 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def dcConv():
        A=Benchmarks.DCConv.A
        B=Benchmarks.DCConv.B
        mode='.'
        E={
        (1,1): [0.9,1.1],
        }
        IS=np.array([
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> DC Conv <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def spaceCraftRndzvs():
        A=Benchmarks.SpaceCraftRndzvs.A
        B=Benchmarks.SpaceCraftRndzvs.B
        mode='.'
        E={
        (2,0): [0.9,1.1],
        (2,3): [0.9,1.1],
        (3,2): [0.9,1.1],
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> SpaceCraft Rndzvs <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def holesCXc():
        A=Benchmarks.HolesCXc.A
        B=Benchmarks.HolesCXc.B
        mode='.'
        E={
        (0,3): [0.9,1.1],
        (1,2): [0.9,1.1],
        (3,2): [0.9,1.1],
        (4,3): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Holes CXc <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def motorTransmission1():
        A=Benchmarks.MotorTransmission1.A
        B=Benchmarks.MotorTransmission1.B
        mode='.'
        E={
        (0,6): [0.9,1.1],
        (1,6): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Motor Transmission 1 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")

    def motorTransmission2():
        A=Benchmarks.MotorTransmission2.A
        B=Benchmarks.MotorTransmission2.B
        mode='.'
        E={
        (4,0): [0.9,1.1],
        (4,1): [0.9,1.1]
        }
        IS=np.array([
        [1],
        [1],
        [1],
        [1],
        [1]
        ])
        T=2000
        mList=['Split','Interval']
        print("\n>>>>>>>>>>> Motor Transmission 2 <<<<<<<<<<<<<\n")
        rs=ReachSet(A,B,mode,E,IS,T,mList)
        rs.compareReachSets()
        print("\n==========================================\n")




#---------------------

Driver.spaceCraftRndzvs()
