import os,sys
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT) # Set the project path

from lib.MoULDySEngine import * # Importing all the functionalities of MoULDyS.


class Aircraft:


    def getDynamics(ohm=1):
        #ohm=1.5
        A=np.array([
        [0,0,1,0],
        [0,0,0,1],
        [0,0,0,-ohm],
        [0,0,ohm,0]
        ])
        B=np.array([])
        h=0.01
        mode='.'

        P=1
        Er={
        (2,3): [1-(P/100),1+(P/100)],
        (3,2): [1-(P/100),1+(P/100)]
        }

        unsafe1=[(-1000,-49.5),(-1000,1000),(-1000,1000),(-1000,1000)]
        unsafe2=[(11,1000),(-1000,1000),(-1000,1000),(-1000,1000)]

        unsafeList=[unsafe1,unsafe2]

        return (A,B,Er,mode,unsafeList,h)

    def getReachSetX1X2():
        '''
        Returns the reachable set for x1 and x2
        '''

        C=[0,0,0,0]
        V=np.array([
        [1,0,0,0],
        [0,1,0,0],
        [0,0,1,0],
        [0,0,0,1]
        ])
        P=[(1.1,1.11),(1.1,1.11),(20,20.1),(20,20.1)]
        initialSet=(C,V,P)
        T=400

        (A,B,Er,mode,unsafeList,h)=Aircraft.getDynamics(ohm=1)
        A=Aircraft.createMatrix(A,B,'.',0.01)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rs=Split(A,Er,initialSet,T)
        (reachORS,reachRS)=rs.getReachableSetAllList()
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        print(len(reachORS))
        #exit(0)

        VisualizeAircraft.vizX1X2(reachRS,reachORS)


    def offlineMonitor():

        (A,B,Er,mode,unsafeList,h)=Aircraft.getDynamics(ohm=1)

        ######### Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tpTS='uncertain'



        ######### Encode Initial Set ########
        initialSetInt=[(1.1,1.11),(1.1,1.11),(20,20.1),(20,20.1)]


        ######### Encode Other Parameters ########
        pr1=5
        pr2=10
        delta1=2
        delta2=10
        th1=0

        T=631 # Time upto which logs are generated

        (log,actualBehavior)=mEngine.genLog(initialSetInt,T,pr1,delta2)

        reachSets=mEngine.offlineMonitor(log)

        mEngine.vizMonitor(reachSets,log,None,T,th1,fname="viz_aircraft",vizCoverage=20)











if True:
    #print("Try")
    #Aircraft.getReachSetX1X2()
    Aircraft.offlineMonitor()
    #Aircraft.onlineMonitor()
