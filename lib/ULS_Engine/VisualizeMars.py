import pickle
import matplotlib.pyplot as plt
import statistics as stat
import numpy as np
import seaborn as sns
import pandas as pd
from matplotlib.patches import Ellipse
import sys,os
import matplotlib
from gurobipy import *
from StarOperations import *

from Parameters import *

class VisualizeRS:
    '''
    Visualize a reachable set
    '''
    def vizRS(appx,th1,th2,fname="reachSet"):
        '''
        Visualize the given reachable set
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        (X,Y)=VizRS.getPlotsLineFine(appx,th1,th2)
        plt.plot(X,Y,'bo')
        if False:
            plt.savefig(ROVER_RESULTS+fname)
            plt.close()
        else:
            plt.show()

    def vizTraj(traj,th1,th2,fname="trajs"):
        '''
        Visualize a trajectory, i.e. a set of reachable sets
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        pltList=VizRS.getPlotsLineFineList(traj,th1,th2)
        for pt in pltList:
            Xt=pt[0] # RS X
            Yt=pt[1] # RS Y
            plt.plot(Xt,Yt)
        if False:
            plt.savefig(ROVER_RESULTS+fname)
            plt.close()
        else:
            plt.show()

    def vizTrajRobotCloud(traj_rbt,traj_cld,th1,th2,fname="trajsRC"):
        '''
        Visualize a trajectory, i.e. a set of reachable sets
        for both robot and cloud
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        pltList_rbt=VizRS.getPlotsLineFineList(traj_rbt,th1,th2)
        pltList_cld=VizRS.getPlotsLineFineList(traj_cld,th1,th2)

        for pt in pltList_cld:
            Xt=pt[0] # RS X
            Yt=pt[1] # RS Y
            plt.plot(Xt,Yt,'g.')

        for pt in pltList_rbt:
            Xt=pt[0] # RS X
            Yt=pt[1] # RS Y
            plt.plot(Xt,Yt,'b.')

        if False:
            plt.savefig(ROVER_RESULTS+fname)
            plt.close()
        else:
            plt.show()

    def vizTrajRobotCloudData(traj_rbt,traj_cld,title,th1=0,th2=1,fname="data_info"):
        '''
        Visualize a trajectory, i.e. a set of reachable sets
        for both robot and cloud
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.xlabel("State "+str(th1))
        plt.ylabel("State "+str(th2))
        plt.title(title)

        pltList_rbt=VizRS.getPlotsLineFineList(traj_rbt,th1,th2)
        pltList_cld=VizRS.getPlotsLineFineList(traj_cld,th1,th2)

        for pt in pltList_cld:
            Xt=pt[0] # RS X
            Yt=pt[1] # RS Y
            plt.plot(Xt,Yt,'g.')

        for pt in pltList_rbt:
            Xt=pt[0] # RS X
            Yt=pt[1] # RS Y
            plt.plot(Xt,Yt,'b.')

        if True:
            plt.savefig(ROVER_RESULTS+fname)
            plt.close()
        else:
            plt.show()

    def vizUnsafeSet(unsafelist,pointClouds):
        n=len(unsafelist[0][0])
        for q in it.combinations(list(range(n)),2):
            th1=q[0]
            th2=q[1]
            #print(th1,th2)
            plt.xlabel("State "+str(th1))
            plt.ylabel("State "+str(th2))
            pltList_unsafe=VizRS.getPlotsLineFineList(unsafelist,th1,th2)

            for pt in pltList_unsafe:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'g.')


            # Plot the scatter points

            for pointCloud in pointClouds:
                X=[p[th1] for p in pointCloud]
                Y=[p[th2] for p in pointCloud]
                plt.scatter(X, Y, c='r',s=100)

            if False:
                plt.savefig(ROVER_RESULTS+fname)
                plt.close()
            else:
                plt.show()

    def vizTrajRobotCloudAllSides(traj_rbt,traj_cld,unsafelist,pointClouds,fname="trajsRC"):
        '''
        Visualize a trajectory, i.e. a set of reachable sets
        for both robot and cloud
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        n=len(traj_cld[0][0])
        for q in it.combinations(list(range(n)),2):
            th1=q[0]
            th2=q[1]
            #print(th1,th2)
            plt.xlabel("State "+str(th1))
            plt.ylabel("State "+str(th2))
            pltList_rbt=VizRS.getPlotsLineFineList(traj_rbt,th1,th2)
            pltList_cld=VizRS.getPlotsLineFineList(traj_cld,th1,th2)
            pltList_unsafe=VizRS.getPlotsLineFineList(unsafelist,th1,th2)

            # Plot the scatter points
            for pointCloud in pointClouds:
                X=[p[th1] for p in pointCloud]
                Y=[p[th2] for p in pointCloud]
                plt.scatter(X, Y, c='b',s=2)

            for pt in pltList_unsafe:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'r.')

            for pt in pltList_cld:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'g.')

            for pt in pltList_rbt:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'b.')

            if False:
                plt.savefig(ROVER_RESULTS+fname)
                plt.close()
            else:
                plt.show()

    def vizTrajRobotCloudExpAllSides(traj_rbt,traj_cld,traj_exp,unsafelist,title,fname="trajsRC_Exp"):
        '''
        Visualize a trajectory, i.e. a set of reachable sets
        for both robot and cloud
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        plt.title(title)
        n=len(traj_cld[0][0])
        for q in it.combinations(list(range(n)),2):
            th1=q[0]
            th2=q[1]
            #print(th1,th2)
            plt.xlabel("State "+str(th1))
            plt.ylabel("State "+str(th2))
            pltList_rbt=VizRS.getPlotsLineFineList(traj_rbt,th1,th2)
            pltList_cld=VizRS.getPlotsLineFineList(traj_cld,th1,th2)
            pltList_expCld=VizRS.getPlotsLineFineList(traj_exp,th1,th2)
            pltList_unsafe=VizRS.getPlotsLineFineList(unsafelist,th1,th2)

            for pt in pltList_unsafe:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'r.')

            for pt in pltList_expCld:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'c.')

            for pt in pltList_cld:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'k.')

            for pt in pltList_rbt:
                Xt=pt[0] # RS X
                Yt=pt[1] # RS Y
                plt.plot(Xt,Yt,'b.')

            if False:
                plt.savefig(ROVER_RESULTS+fname)
                plt.close()
            else:
                plt.show()

    def vizBloating(star,bloatedStar):
        '''
        This is just for debugging purposes
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        n=len(star[0])
        for q in it.combinations(list(range(n)),2):
            th1=q[0]
            th2=q[1]
            plt.xlabel("State "+str(th1))
            plt.ylabel("State "+str(th2))
            (X,Y)=VizRS.getPlotsLineFine(star,th1,th2)
            (X_bloat,Y_bloat)=VizRS.getPlotsLineFine(bloatedStar,th1,th2)
            plt.plot(X,Y,'bo')
            plt.plot(X_bloat,Y_bloat,'co')
            plt.show()

    def vizBloating3(robot,cloud,bloatedRobot):
        '''
        This is just for debugging purposes
        '''
        plt.autoscale(enable=True, axis='both', tight=False)
        n=len(robot[0])
        for q in it.combinations(list(range(n)),2):
            th1=q[0]
            th2=q[1]
            plt.xlabel("State "+str(th1))
            plt.ylabel("State "+str(th2))
            (X_r,Y_r)=VizRS.getPlotsLineFine(robot,th1,th2)
            (X_c,Y_c)=VizRS.getPlotsLineFine(cloud,th1,th2)
            (X_rb,Y_rb)=VizRS.getPlotsLineFine(bloatedRobot,th1,th2)
            plt.plot(X_r,Y_r,'bo')
            plt.plot(X_rb,Y_rb,'co')
            plt.plot(X_c,Y_c,'ro')
            plt.show()

class VisualizeMarsRover:
    '''
    Visualization APIs related to Mars Rover Data
    '''
    def getPlotsLineFine(star,theta1,theta2):
        '''
        Returns the list of points (x,y)
        for the reachable set st1
        '''

        C=star[0]
        V=star[1]
        P=star[2]
        X_list=[]
        Y_list=[]

        sv=V.shape[0]
        aS=V.shape[1]

        semiDefFlag=False
        model = Model("qp")
        model.setParam( 'OutputFlag', False )

        # Create Predicate Variables
        predVars=[]
        for i in range(aS):
            name="Pred"+str(i)
            predVars.append(model.addVar(-GRB.INFINITY,GRB.INFINITY,name=name,vtype='C'))
        #-----------------------

        # Axes Variables
        X=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="X",vtype='C')
        Y=model.addVar(-GRB.INFINITY,GRB.INFINITY,name="Y",vtype='C')
        #-------------------------

        # Create the Star Constraints
        objX=0
        for i in range(aS):
            objX=objX+(predVars[i]*V[theta1][i])
        objX=C[theta1]+objX

        objY=0
        for i in range(aS):
            objY=objY+(predVars[i]*V[theta2][i])
        objY=C[theta2]+objY

        model.addConstr(X==objX,"X Axis")
        model.addConstr(Y==objY,"Y Axis")
        #-----------------------------------

        # Predicate Constraints
        for i in range(aS):
            a=P[i][0]
            b=P[i][1]
            if a==b:
                model.addConstr(predVars[i]==a,name)
            else:
                model.addConstr(predVars[i]>=min(a,b),name+".1")
                model.addConstr(predVars[i]<=max(a,b),name+".2")
        #-----------------------------------


        # Quadrant Specific Constraints

        # 1st Quadrant

        obj=X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(900):
            an=a/10
            if an==90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*(X),"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        #'''
        # 2nd Quadrant

        obj=X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(0,-900,-1):
            an=a/10
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 3rd Quadrant

        obj=-X-Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(-900,-1800,-1):
            an=a/10
            if an==-90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------

        # 4th Quadrant

        obj=-X+Y
        model.setObjective(obj,GRB.MAXIMIZE)

        for a in range(900,1800):
            an=a/10
            if an==90:
                model.addConstr(X==0,"Angle")
            else:
                m=math.tan(math.radians(an))
                model.addConstr(Y==m*X,"Angle")
            try:
                model.optimize()
                #model.write("dump.bas")
                status = model.Status
                if status==GRB.Status.UNBOUNDED:
                    print("UNBOUNDED ")
                else:
                    if status == GRB.Status.INF_OR_UNBD or \
                       status == GRB.Status.INFEASIBLE  or \
                       status == GRB.Status.UNBOUNDED:
                        0;
                        #print('**The model cannot be solved because it is infeasible or unbounded**')
                    else:
                        xVal=model.getVarByName("X").x
                        yVal=model.getVarByName("Y").x
                        X_list.append(xVal)
                        Y_list.append(yVal)
            except:
                semiDefFlag=True

            if semiDefFlag==True:
                print("Shoot!!")

            semiDefFlag=False
            model.remove(model.getConstrByName("Angle"))
        #-----------------------------



        #------------------------------

        #print(X_list,Y_list)
        #exit(0)
        #'''

        return (X_list,Y_list)

    def vizMarsPointClouds(obs, terrain_img, fname="PointCloud", xlim=(0,224), ylim=(0,224), title_str="Point Clouds"):
        '''
        Modified from: https://bitbucket.org/nakanoya/tasknet-icra2021/src/master/Motion_Planning/utils/visualize_utils.py
        '''
        xmin = xlim[0]
        xmax = xlim[1]

        ymin = ylim[0]
        ymax = ylim[1]

        fig = plt.figure()
        # point cloud obstacle
        plt.scatter(obs[:,0], obs[:,1], c='red', s=2)

        linestyle = 'solid'
        color = 'b'

        ax = fig.gca()
        ax.set_xticks(np.arange(xmin,xmax,50))
        ax.set_yticks(np.arange(ymin,ymax,50))
        ax.set_aspect('equal')
        ax.imshow(terrain_img, cmap='gray')
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.grid()
        plt.title(title_str)

        if True:
            plt.savefig(ROVER_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def vizMarsTraj(obs, terrain_img, planned_path_xy, fname="Traj", xlim=(0,224), ylim=(0,224), title_str="Point Clouds"):
        '''
        Modified from: https://bitbucket.org/nakanoya/tasknet-icra2021/src/master/Motion_Planning/utils/visualize_utils.py
        '''
        xmin = xlim[0]
        xmax = xlim[1]

        ymin = ylim[0]
        ymax = ylim[1]

        fig = plt.figure()
        # point cloud obstacle
        plt.scatter(obs[:,0], obs[:,1], c='red', s=2)

        linestyle = 'solid'
        color = 'b'

        for path_xy in planned_path_xy:
            path_x = path_xy[0]
            path_y = path_xy[1]

            plt.plot(path_x, path_y, c=color, marker='o', linewidth=0.5, markersize=0.5, linestyle=linestyle, alpha=1.0)

        # start
        #print(planned_path_xy[0][0], planned_path_xy[1][0])
        plt.plot(planned_path_xy[0][0], planned_path_xy[0][1], c='g', marker='o', markersize=6)
        # end
        plt.plot(planned_path_xy[-1][0], planned_path_xy[-1][1], c='y', marker='s', markersize=6)

        ax = fig.gca()
        ax.set_xticks(np.arange(xmin,xmax,50))
        ax.set_yticks(np.arange(ymin,ymax,50))
        ax.set_aspect('equal')
        ax.imshow(terrain_img, cmap='gray')
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.grid()
        plt.title(title_str)

        if True:
            plt.savefig(ROVER_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def vizMarsTrajMPC(obs, terrain_img, path_ref, spline_path, path_mpc, fname="TrajMPC", xlim=(0,224), ylim=(0,224), title_str="Point Clouds"):
        '''
        Modified from: https://bitbucket.org/nakanoya/tasknet-icra2021/src/master/Motion_Planning/utils/visualize_utils.py
        '''
        xmin = xlim[0]
        xmax = xlim[1]

        ymin = ylim[0]
        ymax = ylim[1]

        fig = plt.figure()
        # point cloud obstacle
        plt.scatter(obs[:,0], obs[:,1], c='red', s=2)

        linestyle = 'solid'
        color = 'b'

        plt.plot(spline_path[0], spline_path[1], c=color, marker='o', linewidth=1, markersize=1, linestyle=linestyle, alpha=1.0)

        for path in path_ref:
            x_ref=[p[0] for p in path]
            y_ref=[p[1] for p in path]
            plt.scatter(x_ref, y_ref, marker='*', c='blue', s=70)

        # start
        plt.plot(path_ref[0][0][0], path_ref[0][0][1], c='g', marker='o', markersize=7)
        # end
        plt.plot(path_ref[-1][-1][0], path_ref[-1][-1][1], c='y', marker='s', markersize=7)

        for path_xy in path_mpc:
            path_x = path_xy[0]
            path_y = path_xy[1]

            plt.plot(path_x, path_y, c='m', marker='o', linewidth=4, markersize=4, linestyle=linestyle, alpha=0.7)


        ax = fig.gca()
        ax.set_xticks(np.arange(xmin,xmax,50))
        ax.set_yticks(np.arange(ymin,ymax,50))
        ax.set_aspect('equal')
        ax.imshow(terrain_img, cmap='gray')
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.grid()
        plt.title(title_str)

        if True:
            plt.savefig(ROVER_RESULTS+"/"+fname, dpi=100, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def vizMarsTrajRS(obs, terrain_img, spline_path, mpc_path, ref_path, traj, fname="TrajRS", xlim=(0,224), ylim=(0,224), title_str="Point Clouds"):
        '''
        Modified from: https://bitbucket.org/nakanoya/tasknet-icra2021/src/master/Motion_Planning/utils/visualize_utils.py
        '''
        xmin = xlim[0]
        xmax = xlim[1]

        ymin = ylim[0]
        ymax = ylim[1]

        fig = plt.figure()

        linestyle = 'solid'
        color = 'b'

        th1=0
        th2=1

        plt.plot(spline_path[0], spline_path[1], c='black', marker='o', linewidth=1, markersize=1, linestyle='-', alpha=1.0)

        for rsList in traj:
            for rs in rsList:
                (X,Y)=VisualizeMarsRover.getPlotsLineFine(rs,th1,th2)
                fgSafe=StarOp.checkIntersectionPoints(rs,obs)
                if fgSafe==True:
                    plt.plot(X,Y,'r.',alpha=0.05)
                else:
                    plt.plot(X,Y,'b.',alpha=0.05)

        for path_xy in mpc_path:
            path_x = path_xy[0]
            path_y = path_xy[1]

            plt.plot(path_x, path_y, c='m', marker='o', linewidth=0.5, markersize=1, linestyle=linestyle, alpha=1.0)

        for path in ref_path:
            x_ref=[p[0] for p in path]
            y_ref=[p[1] for p in path]
            plt.scatter(x_ref, y_ref, marker='*', c='m', s=70, alpha=1)

        # start
        plt.plot(ref_path[0][0][0], ref_path[0][0][1], c='g', marker='o', markersize=7)
        # end
        plt.plot(ref_path[-1][-1][0], ref_path[-1][-1][1], c='y', marker='s', markersize=7)

        # point cloud obstacle
        plt.scatter(obs[:,0], obs[:,1], c='red', s=2)

        ax = fig.gca()
        ax.set_xticks(np.arange(xmin,xmax,50))
        ax.set_yticks(np.arange(ymin,ymax,50))
        ax.set_aspect('equal')
        ax.imshow(terrain_img, cmap='gray')
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.grid()
        plt.title(title_str)

        if True:
            plt.savefig(ROVER_RESULTS+fname, dpi=100, bbox_inches='tight')
            plt.close()
        else:
            plt.show()

    def vizMarsTrajRSComp(obs, terrain_img, spline_path, mpc_path, ref_path, trajTop, trajBot, fname="TrajRSComp", xlim=(0,224), ylim=(0,224), title_str="Point Clouds"):
        '''
        Modified from: https://bitbucket.org/nakanoya/tasknet-icra2021/src/master/Motion_Planning/utils/visualize_utils.py
        '''
        xmin = xlim[0]
        xmax = xlim[1]

        ymin = ylim[0]
        ymax = ylim[1]

        fig = plt.figure()

        linestyle = 'solid'
        color = 'b'

        th1=0
        th2=1

        plt.plot(spline_path[0], spline_path[1], c='black', marker='o', linewidth=1, markersize=1, linestyle='-', alpha=1.0)

        for rsList in trajTop:
            for rs in rsList:
                (X,Y)=VisualizeMarsRover.getPlotsLineFine(rs,th1,th2)
                fgSafe=StarOp.checkIntersectionPoints(rs,obs)
                if fgSafe==True:
                    plt.plot(X,Y,'r.',alpha=0.05)
                else:
                    plt.plot(X,Y,'b.',alpha=0.05)

        for rsList in trajBot:
            for rs in rsList:
                (X,Y)=VisualizeMarsRover.getPlotsLineFine(rs,th1,th2)
                fgSafe=StarOp.checkIntersectionPoints(rs,obs)
                if fgSafe==True:
                    plt.plot(X,Y,'m.',alpha=0.1)
                else:
                    plt.plot(X,Y,'g.',alpha=0.1)

        # point cloud obstacle
        plt.scatter(obs[:,0], obs[:,1], c='red', s=2)

        for path_xy in mpc_path:
            path_x = path_xy[0]
            path_y = path_xy[1]

            plt.plot(path_x, path_y, c='m', marker='o', linewidth=0.5, markersize=1, linestyle=linestyle, alpha=1.0)

        for path in ref_path:
            x_ref=[p[0] for p in path]
            y_ref=[p[1] for p in path]
            plt.scatter(x_ref, y_ref, marker='*', c='m', s=70, alpha=1)

        # start
        plt.plot(ref_path[0][0][0], ref_path[0][0][1], c='g', marker='o', markersize=7)
        # end
        plt.plot(ref_path[-1][-1][0], ref_path[-1][-1][1], c='y', marker='s', markersize=7)

        ax = fig.gca()
        ax.set_xticks(np.arange(xmin,xmax,50))
        ax.set_yticks(np.arange(ymin,ymax,50))
        ax.set_aspect('equal')
        ax.imshow(terrain_img, cmap='gray')
        plt.xlim(xmin, xmax)
        plt.ylim(ymax, ymin)
        plt.grid()
        plt.title(title_str)

        if True:
            plt.savefig(ROVER_RESULTS+fname, dpi=100, bbox_inches='tight')
            plt.close()
        else:
            plt.show()
