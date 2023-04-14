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
            I=np.identity(n1,dtype=np.float)
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

        VisualizeAircraft.vizRS(reachRS,[])


    def offlineCaseStudy(loadResultsFromPickle=True,fig=1):

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

        #pr=pr1
        #delta=delta1
        if loadResultsFromPickle==True:
            if fig==1:
                pickleFileName='aircraftPickle_'+str(pr1)+'pr_'+str(delta1)+'utime'
            elif fig==2:
                pickleFileName='aircraftPickle_'+str(pr1)+'pr_'+str(delta2)+'utime'
            elif fig==3:
                pickleFileName='aircraftPickle_'+str(pr2)+'pr_'+str(delta1)+'utime'
            elif fig==4:
                pickleFileName='aircraftPickle_'+str(pr2)+'pr_'+str(delta2)+'utime'
            else:
                print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid choice!{bcolors.ENDC}")
                exit()
            objects = []
            #pickleFileName='aircraftPickle_'+str(pr)+'pr_'+str(delta)+'utime'
            with open(DATA_PATH+'/aircraft/pickles/'+pickleFileName, 'rb') as handle:
                objects = pickle.load(handle)
            reachSets=(objects['reachSets'],objects['reachSetsTimes'])
            log=Logs(objects['logs'],'uncertain','zono')
        else:
            ############ Generate Logs ############
            (log,actualBehavior)=mEngine.genLog(initialSetInt,T,pr,delta)
            reachSets=mEngine.offlineMonitor(log)

        mEngine.vizMonitor(reachSets,log,None,T,th1,fname="viz_aircraft",vizCoverage=20)


if len(sys.argv)<=1:
    print(f"{bcolors.BOLD}{bcolors.FAIL}No arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)>3:
    print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)==3:
    if sys.argv[1]=='-offline':
        Aircraft.offlineCaseStudy(True,int(sys.argv[2]))
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
elif len(sys.argv)==2:
    if sys.argv[1]=='-behavior':
        Aircraft.getReachSetX1X2()
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
else:
    print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
    exit()
