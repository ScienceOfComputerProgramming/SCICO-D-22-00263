import os,sys
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT) # Set the project path

from lib.MoULDySEngine import * # Importing all the functionalities of MoULDyS.


class Aircraft:

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
        T=631

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

        Visualize.vizRS(reachRS,[])


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


    def onlineMonitor():
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

        (log,actualBehavior)=mEngine.genLog(initialSetInt,T,pr1,0)

        (reachSets,logs)=mEngine.onlineMonitor(actualBehavior)

        mEngine.vizMonitor(reachSets,logs,None,T,th1,fname="viz_aircraft",vizCoverage=20)


    def compareMonitoring():

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

        (logOffline,actualBehavior)=mEngine.genLog(initialSetInt,T,pr1,delta2)

        reachSetsOffline=mEngine.offlineMonitor(logOffline)

        (reachSetsOnline,logsOnline)=mEngine.onlineMonitor(actualBehavior)

        mEngine.vizCompMonitor(reachSetsOffline,logOffline,reachSetsOnline,logsOnline,None,T,th1,fname="viz_aircraft",vizCoverage=VIZ_PER_COVERAGE,tpTS=tpTS)


    def offlineMonitorFile():

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

        pr=pr1
        delta=delta2

        T=631 # Time upto which logs are generated

        loglogFname="aircraft/logs/_5pr_10utime_interval"

        #fname="aircraft/logs/"

        #(log,actualBehavior)=mEngine.genLogFile(initialSetInt,T,fname,"zonotope",pr,delta)

        reachSets=mEngine.offlineMonitorLogFile(loglogFname,tpRep='interval',tpTS='uncertain')

        mEngine.vizMonitorLogFile(reachSets,loglogFname,"interval",T,th1,fname="viz_aircraft",vizCoverage=20,tpTS='uncertain')

    def onlineMonitorFile():

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

        pr=pr1
        delta=delta2

        T=631 # Time upto which logs are generated

        loglogFname="aircraft/logs/_"+str(pr)+"pr_"+str(delta)+"utime_interval"

        #fname="aircraft/logs/"

        #(log,actualBehavior)=mEngine.genLogFile(initialSetInt,T,fname,"zonotope",pr,delta)

        reachSets,logs=mEngine.onlineMonitorBehFile(loglogFname,tp='interval')

        mEngine.vizMonitor(reachSets,logs,None,T,th1,fname="viz_aircraft",vizCoverage=20)

    def compareMonitorFile():

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

        pr=pr1
        delta=delta2

        T=631 # Time upto which logs are generated

        loglogFname="aircraft/logs/_"+str(pr)+"pr_"+str(delta)+"utime_interval"

        reachSetsOffline=mEngine.offlineMonitorLogFile(loglogFname,tpRep='interval',tpTS='uncertain')

        #fname="aircraft/logs/"

        #(log,actualBehavior)=mEngine.genLogFile(initialSetInt,T,fname,"zonotope",pr,delta)

        reachSetsOnline,logsOnline=mEngine.onlineMonitorBehFile(loglogFname,tp='interval')

        mEngine.vizCompMonitorLogFile(reachSetsOffline,loglogFname,reachSetsOnline,logsOnline,"interval",T,th1,fname="viz_test",vizCoverage=VIZ_PER_COVERAGE,tpTS='uncertain')










if True:
    #print("Try")
    #Aircraft.getReachSetX1X2()
    #Aircraft.offlineMonitor()
    #Aircraft.onlineMonitor()
    #Aircraft.compareMonitoring()
    #Aircraft.offlineMonitorFile()
    #Aircraft.onlineMonitorFile()
    Aircraft.compareMonitorFile()
