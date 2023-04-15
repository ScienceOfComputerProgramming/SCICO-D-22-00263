import os,sys
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT) # Set the project path

from lib.MoULDySEngine import * # Importing all the functionalities of MoULDyS.

class Anesthesia:

    def getDynamics():
        ######### Step 1: Define Dynamics Matrices ########
        #With weight of 25 Kg
        weight=25
        v1=458.4*weight
        k10=0.1527*pow(weight,-0.3)
        k12=0.114
        k13=0.0419
        k21=0.055
        k31=0.0033
        kd=40 #with Td=20s
        A=np.array([
        [-(k10+k12+k13),k12,k13,0],
        [k21,-k21,0,0],
        [k31,0,-k31,0],
        [kd,0,0,-kd]
        ])
        B=np.array([
        [1/v1],
        [0],
        [0],
        [0]
        ])

        ######### Step 2: Mode ########
        mode='.' # Continuous mode
        h=0.01 # Discretization parameter

        ######### Step 3: Encode Uncertainties ########
        p=0.8 # 0.8
        Er={
        (0,0): [1-(p/2000),1+(p/2000)],
        (0,4): [100/(100+p),100/(100-p)]
        }

        ######### Step 4: Unsafe Set ########
        '''
        The behavior is marked as unsafe if: c_p<=1 OR c_p>=6. This encoded by `unsafe`.
        '''
        #unsafe1=[(-np.inf,1),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)]
        #unsafe2=[(6,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)]
        unsafe1=[(-1000,1),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        unsafe2=[(6,1000),(-1000,1000),(-1000,1000),(-1000,1000),(-1000,1000)]
        unsafeList=[unsafe1,unsafe2]

        return (A,B,Er,mode,unsafeList,h)

    def offlineCaseStudy(fig=1):

        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=Anesthesia.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'


        '''
        ######### Encode Initial Set (optional) ########
        initialSet1=[(3,4),(3,4),(4,5),(3,4),(2,5)] # Initial Set, represented as an interval
        initialSet2=[(2,4),(3,6),(3,6),(2,4),(2,10)] # Initial Set, represented as an interval
        pr1=20
        pr2=40
        pr3=5
        T=2000 # Time upto which logs are generated

        tp='interval' # Use tp='zonotope' for creating logs with zonotopes
        if True:
            fname='anesthesia/anesthesia_recreate_ISlarge'
            T=2000 # Time upto which logs are generated
            np.random.seed()
            (log,actualBehavior)=mEngine.genLogFile(initialSet2,T,fname,tp,pr2)
        '''


        ######### Step 6: Log File ########
        if fig==1:
            logFname='anesthesia_recreate_ISsmall_20_'+tp
        elif fig==2:
            logFname='anesthesia_recreate_ISlarge_20_'+tp
        elif fig==3:
            logFname='anesthesia_recreate_ISsmall_40_'+tp
        elif fig==4:
            logFname='anesthesia_recreate_ISlarge_40_'+tp
        else:
            print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid choice!{bcolors.ENDC}")
            exit()

        '''
        if fig==1:
            (logs,ab)=mEngine.genLog(initialSet1,T,pr1)
        elif fig==2:
            (logs,ab)=mEngine.genLog(initialSet1,T,pr2)
        elif fig==3:
            (logs,ab)=mEngine.genLog(initialSet2,T,pr1)
        elif fig==4:
            (logs,ab)=mEngine.genLog(initialSet2,T,pr2)
        '''


        ######### Step 7: Perform Offline Monitoring ########
        reachSets=mEngine.offlineMonitorLogFile("anesthesia/"+logFname,tp)

        ######### Step 7: Visualize the results of monitoring ########
        T=2000 # Time step upto which we want to visualize
        th1=0 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizMonitorLogFile(reachSets,"/anesthesia/"+logFname,tp,T,th1,"anesthesia_offline_monitor_fig"+str(fig),vizCoverage=vizCov)

    def onlineCaseStudy():
        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=Anesthesia.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'

        ######### Step 6: Behavior File ########
        logFname='anesthesia_recreate_ISsmall_40_interval'

        ######### Step 7: Perform Online Monitoring ########
        (reachSets,logs)=mEngine.onlineMonitorBehFile("/anesthesia/"+logFname,tp)

        T=2000 # Time step upto which we want to visualize
        th1=0 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizMonitor(reachSets,logs,tp,T,th1,"anesthesia_online_monitor",vizCoverage=vizCov)

    def compareCaseStudy():
        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=Anesthesia.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'

        ######### Step 6: Log and Behavior File ########
        logFname='anesthesia_recreate_ISsmall_5_interval'

        ######### Step 7: Compare Offline and Online Monitoring ########
        reachSets=mEngine.offlineMonitorLogFile("anesthesia/"+logFname,tp)
        (reachSetsOnline,logsOnline)=mEngine.onlineMonitorBehFile("/anesthesia/"+logFname,tp)

        T=2000 # Time step upto which we want to visualize
        th1=0 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizCompMonitorLogFile(reachSets,"/anesthesia/"+logFname,reachSetsOnline,logsOnline,tp,T,th1,"anesthesia_comp_monitor",vizCov)




if len(sys.argv)<=1:
    print(f"{bcolors.BOLD}{bcolors.FAIL}No arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)>3:
    print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)==3:
    if sys.argv[1]=='-offline':
        Anesthesia.offlineCaseStudy(int(sys.argv[2]))
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
elif len(sys.argv)==2:
    if sys.argv[1]=='-compare':
        Anesthesia.compareCaseStudy()
    elif sys.argv[1]=='-online':
        Anesthesia.onlineCaseStudy()
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
