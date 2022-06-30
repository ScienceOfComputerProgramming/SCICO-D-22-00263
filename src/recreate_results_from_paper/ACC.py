import os,sys
PROJECT_ROOT = os.environ['MNTR_ROOT_DIR']
sys.path.append(PROJECT_ROOT) # Set the project path

from lib.MoULDySEngine import * # Importing all the functionalities of MoULDyS.

class ACC:

    def getDynamics():
        ######### Step 1: Define Dynamics Matrices ########
        m=1370
        g=9.8
        f1=1.2567
        F=-1.2156
        aL=0.5
        A=np.array([
        [f1/m,0,0,F],
        [-1,0,1,0],
        [0,0,0,aL],
        [0,0,0,0]
        ])
        B=np.array([])

        ######### Step 2: Mode ########
        mode='.' # Continuous mode
        h=0.01 # Discretization parameter

        ######### Step 3: Encode Uncertainties ########
        Er={
        (0,3): [-0.6,2.46],
        (2,3): [-0.9,0.6]
        }

        ######### Step 4: Unsafe Set ########
        '''
        The behavior is marked as unsafe if: h<=0.5. This encoded by `unsafe`.
        '''
        unsafe1=[(-np.inf,np.inf),(-np.inf,0.5),(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)]
        unsafeList=[unsafe1]

        return (A,B,Er,mode,unsafeList,h)

    def offlineCaseStudy(fig=1):

        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=ACC.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'


        '''
        ######### Encode Initial Set (optional) ########
        initialSet1=[(15,15.01),(3,3.03),(14.9,15),(1,1)] # Initial Set, represented as an interval
        initialSet2=[(15,15.1),(3,3.5),(14.9,15.1),(1,1)] # Initial Set, represented as an interval
        pr1=20
        pr2=40
        pr3=14
        T=2000 # Time upto which logs are generated

        tp='interval' # Use tp='zonotope' for creating logs with zonotopes
        if True:
            fname='acc/acc_recreate_ISsmall'
            T=2000 # Time upto which logs are generated
            (log,actualBehavior)=mEngine.genLogFile(initialSet1,T,fname,tp,pr3)
            exit()
        '''


        ######### Step 6: Log File ########
        if fig==1:
            logFname='acc_recreate_ISsmall_20_'+tp
        elif fig==2:
            logFname='acc_recreate_ISlarge_20_'+tp
        elif fig==3:
            logFname='acc_recreate_ISsmall_40_'+tp
        elif fig==4:
            logFname='acc_recreate_ISlarge_40_'+tp
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
        reachSets=mEngine.offlineMonitorLogFile("acc/"+logFname,tp)

        ######### Step 7: Visualize the results of monitoring ########
        T=2000 # Time step upto which we want to visualize
        th1=1 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizMonitorLogFile(reachSets,"/acc/"+logFname,tp,T,th1,"acc_offline_monitor",vizCoverage=vizCov)

    def onlineCaseStudy():
        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=ACC.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'

        ######### Step 6: Behavior File ########
        logFname='acc_recreate_ISsmall_40_interval'

        ######### Step 7: Perform Online Monitoring ########
        (reachSets,logs)=mEngine.onlineMonitorBehFile("/acc/"+logFname,tp)

        T=2000 # Time step upto which we want to visualize
        th1=1 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizMonitor(reachSets,logs,tp,T,th1,"/acc/"+logFname,vizCoverage=vizCov)

    def compareCaseStudy():
        ######### Step 1-4 ########
        (A,B,Er,mode,unsafeList,h)=ACC.getDynamics()

        ######### Step 5: Instantiate the MoULDyS engine ########
        mEngine=MoULDyS(A,B,Er,mode,unsafeList,h)
        tp='interval'

        ######### Step 6: Log and Behavior File ########
        logFname='acc_recreate_ISsmall_14_interval'

        ######### Step 7: Compare Offline and Online Monitoring ########
        reachSets=mEngine.offlineMonitorLogFile("acc/"+logFname,tp)
        (reachSetsOnline,logsOnline)=mEngine.onlineMonitorBehFile("/acc/"+logFname,tp)

        T=2000 # Time step upto which we want to visualize
        th1=0 # State variable that is to be visualized
        vizCov=20 # Percentage of reachable sets to be visualized. Note: Visualizing all reachable sets is expensive.
        #Note: Visualization takes time!
        mEngine.vizCompMonitorLogFile(reachSets,"/acc/"+logFname,reachSetsOnline,logsOnline,tp,T,th1,"viz_test",vizCov)




if len(sys.argv)<=1:
    print(f"{bcolors.BOLD}{bcolors.FAIL}No arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)>3:
    print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
    exit()
elif len(sys.argv)==3:
    if sys.argv[1]=='-offline':
        ACC.offlineCaseStudy(int(sys.argv[2]))
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
elif len(sys.argv)==2:
    if sys.argv[1]=='-compare':
        ACC.compareCaseStudy()
    elif sys.argv[1]=='-online':
        ACC.onlineCaseStudy()
    else:
        print(f"{bcolors.BOLD}{bcolors.FAIL}Invalid arguments provided!{bcolors.ENDC}")
        exit()
