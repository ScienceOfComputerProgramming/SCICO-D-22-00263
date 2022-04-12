'''
APIs for Offloading with Mars Rover
'''
from VisualizeMars import *
#from RobotCloudModel import *
from PythonRobotics.PathTracking.model_predictive_speed_and_steer_control.model_predictive_speed_and_steer_control import *
import sys,os,copy
import numpy as np
import time
from Parameters import *
from BicycleModel import *

class Rover:
    '''
    APIs for Offloading on the Mars Rover data
    '''

    def __init__(self,env_id=ENV_ID,mars_point_cloud_fname=MARS_POINT_CLOUD):
        self.mars_point_cloud_fname=mars_point_cloud_fname
        self.env_id=env_id
        self.image=[]
        self.obs=[]
        self.loadData()

    def loadData(self):
        '''
        Load the data in objects from the file
        '''
        dataset = np.load(self.mars_point_cloud_fname)
        self.image  = dataset['image'][self.env_id]
        self.obs = dataset['point_cloud'][self.env_id]


    def getReachSets(self):
        '''
        Get Reachable Sets of Robot and Cloud
        '''
        VisualizeMarsRover.vizMarsPointClouds(self.obs,self.image)

        if self.env_id==2:
            path=[
            [(75,25),(52,70),(55,110),
            (100,100),
            (70,135),(80,190),(140,200),
            (136,160),
            (150,100),(100,50),(85,40)]
            ]
        elif self.env_id==0:
            path=[
            [(25,50),(85,70),(75,118),(100,142),(71,175),(15,167)]
            ]
        elif self.env_id==1:
            case=5
            if case==2:
                path=[
                [(20,115),(50,125),(100,129),(112,160),(165,170)]
                ]
            elif case==3:
                path=[
                [(20,115),(65,79),(100,85),(125,128),(165,170)]
                ]
            elif case==4:
                path=[
                [(20,115),(60,60),(128,76),(165,170)]
                ]
            elif case==5:
                path=[
                [(20,115),(60,80),(100,80),(120,110),(120,140),(40,115)]
                ]
        else:
            print("Invalid Environment Chosen!")
            exit(0)
        #VisualizeMarsRover.vizMarsTraj(self.obs,self.image,path)
        (T, X, Y, yawList, V, D, A, spline_path) = doMPC(path)

        path_mpc=[]
        for x,y in zip(X,Y):
            path_mpc.append((x,y))
        VisualizeMarsRover.vizMarsTrajMPC(self.obs,self.image,path,spline_path,path_mpc)

        print(">> STATUS: Computing Reachable Sets . . .")
        time_taken=time.time()
        rsList=BicycleModel.getReachSets(X, Y, yawList, V, D, A)
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Reachable Sets Computed!")
        #exit(0)
        VisualizeMarsRover.vizMarsTrajRS(self.obs,self.image,spline_path,path_mpc,path,rsList)

    def getReachSetsTopBot(self):
        '''
        Get Reachable Sets of Robot and Cloud, comparing the top and bottom.
        '''
        VisualizeMarsRover.vizMarsPointClouds(self.obs,self.image)

        if self.env_id==2:
            path=[
            [(75,25),(52,70),(55,110),
            (100,100),
            (70,135),(80,190),(140,200),
            (136,160),
            (150,100),(100,50),(85,40)]
            ]
        elif self.env_id==0:
            path=[
            [(25,50),(85,70),(75,118),(100,142),(71,175),(15,167)]
            ]
        elif self.env_id==1:
            case=5
            if case==2:
                path=[
                [(20,115),(50,125),(100,129),(112,160),(165,170)]
                ]
            elif case==3:
                path=[
                [(20,115),(65,79),(100,85),(125,128),(165,170)]
                ]
            elif case==4:
                path=[
                [(20,115),(60,60),(128,76),(165,170)]
                ]
            elif case==5:
                path=[
                [(20,115),(60,80),(100,80),(120,110),(120,140),(40,115)]
                ]
        else:
            print("Invalid Environment Chosen!")
            exit(0)
        #VisualizeMarsRover.vizMarsTraj(self.obs,self.image,path)
        (T, X, Y, yawList, V, D, A, spline_path) = doMPC(path)

        path_mpc=[]
        for x,y in zip(X,Y):
            path_mpc.append((x,y))
        VisualizeMarsRover.vizMarsTrajMPC(self.obs,self.image,path,spline_path,path_mpc)

        (Er_top,Er_bot)=self.getErrors()
        rsListTop=BicycleModel.getReachSetsEr(X, Y, yawList, V, D, A, Er_top)
        rsListBot=BicycleModel.getReachSetsEr(X, Y, yawList, V, D, A, Er_bot)

        print(">> STATUS: Reachable Sets Computed!")

        VisualizeMarsRover.vizMarsTrajRSComp(self.obs,self.image,spline_path,path_mpc,path,rsListTop,rsListBot)


    def getCellOrder(self):
        '''
        Get cell orders
        '''
        v=0.1
        phi=0.1
        delta=0.1
        print(">> STATUS: Computing Cell Order . . .")
        time_taken=time.time()
        cellOrder=BicycleModel.getCellOrderBasedError(v,phi,delta)
        time_taken=time.time()-time_taken
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Cell Order Computed!")
        return cellOrder

    def getErrors(self):
        '''
        Based on the ouput we got from `self.getCellOrder()`,
        some perturbation of the matrix were picked.
        '''
        p=10
        Er_top_0={
        (0,4): [1-(p/100),1+(p/100)],
        (0,5): [1-(p/100),1+(p/100)],
        (1,2): [1-(p/100),1+(p/100)],
        (1,3): [1-(p/100),1+(p/100)]
        }

        Er_bot_0={
        (0,3): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (1,6): [1-(p/100),1+(p/100)],
        (1,5): [1-(p/100),1+(p/100)]
        }

        Er_top={
        (0,6): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (1,2): [1-(p/100),1+(p/100)],
        (1,3): [1-(p/100),1+(p/100)]
        }

        Er_bot={
        (0,3): [1-(p/100),1+(p/100)],
        (0,2): [1-(p/100),1+(p/100)],
        (1,6): [1-(p/100),1+(p/100)],
        (1,5): [1-(p/100),1+(p/100)]
        }

        return Er_top,Er_bot

    def getRobustnessMetric(self):
        '''
        Returns the Robustness Metric, amount of error
        it can tolerate
        '''
        p=9
        if self.env_id==2:
            path=[
            [(75,25),(52,70),(55,110),
            (100,100),
            (70,135),(80,190),(140,200),
            (136,160),
            (150,100),(100,50),(85,40)]
            ]
        elif self.env_id==0:
            path=[
            [(25,50),(85,70),(75,118),(100,142),(71,175),(15,167)]
            ]
        elif self.env_id==1:
            case=5
            if case==2:
                path=[
                [(20,115),(50,125),(100,129),(112,160),(165,170)]
                ]
            elif case==3:
                path=[
                [(20,115),(65,79),(100,85),(125,128),(165,170)]
                ]
            elif case==4:
                path=[
                [(20,115),(60,60),(128,76),(165,170)]
                ]
            elif case==5:
                path=[
                [(20,115),(60,80),(100,80),(120,110),(120,140),(40,115)]
                ]
        else:
            print("Invalid Environment Chosen!")
            exit(0)
        #VisualizeMarsRover.vizMarsTraj(self.obs,self.image,path)
        (T, X, Y, yawList, V, D, A, spline_path) = doMPC(path)

        path_mpc=[]
        for x,y in zip(X,Y):
            path_mpc.append((x,y))

        step=0.5
        print(">> STATUS: Computing Robustness Metric . . .")
        time_taken=time.time()
        rsListOld=[]
        while p<=11:
            Er={
            #(0,0): [1-(p/100),1+(p/100)],
            (0,2): [1-(p/100),1+(p/100)],
            (0,3): [1-(p/100),1+(p/100)],
            (0,6): [1-(p/100),1+(p/100)],
            #(1,1): [1-(p/100),1+(p/100)],
            (1,2): [1-(p/100),1+(p/100)],
            (1,3): [1-(p/100),1+(p/100)],
            (1,6): [1-(p/100),1+(p/100)],
            }
            rsList=BicycleModel.getReachSetsEr(X, Y, yawList, V, D, A, Er)
            safe=BicycleModel.isSafe(rsList,self.obs)
            if safe==False:
                p=p-step
                break;
            p=p+step
            rsListOld=copy.copy(rsList)
        time_taken=time.time()-time_taken
        print("\tRobustness Metric: ",p)
        print("\tTime Taken: ",time_taken)
        print(">> STATUS: Robustness Metric Computed")
        print(">> STATUS: Visualizing Safe Reachable Sets  . . .")
        VisualizeMarsRover.vizMarsTrajRS(self.obs,self.image,spline_path,path_mpc,path,rsListOld,fname="RobMetTrajRS")
        print(">> STATUS: Safe Reachable Sets Visualized!")



if True:
    rov=Rover()
    #rov.getReachSetsTopBot()
    #rov.getReachSets()
    #rov.getCellOrder()
    rov.getRobustnessMetric()
