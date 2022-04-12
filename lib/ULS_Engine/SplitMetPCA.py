'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Linear Dynamical System: dot{x} = (A+E)x; where E is the perturbation.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import time
import numpy.linalg as LA

from VisualizationReachSet import *
from ComputeU import *
from SamplingPCA import *

INTERVAL=300


class SplitPCA:
    '''
    Experimental approach to select basis vectors based on PCA
    '''

    def __init__(self,A,E,theta,T):
        self.A=A # Linear system (without perturbation)
        self.Er=E
        '''
        Following dictionary data-structure has been used to represent
        a error matrix:
        {
            (i,j): [a,b]
        }
        Indicating array[i][j] has a perturbation of -(1-a)% to (b-1)%
        For eg: [0.9,1,1] means a perturbation of -10% and +10%
        '''
        self.Theta=theta # The initial set
        self.T=T #Max number of steps
        self.n=A.shape[0]
        self.Ac=self.computeCenter()

    def computeCenter(self):
        '''
        Computes the center point matrix of the interval matrix self.A
        '''
        Ac=np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                if (i,j) in self.Er:
                    a=float(self.A[i][j]*self.Er[(i,j)][0])
                    b=float(self.A[i][j]*self.Er[(i,j)][1])
                    c=(a+b)/2
                    Ac[i][j]=c
                else:
                    Ac[i][j]=self.A[i][j]
        return Ac

    def printReachableSet(self,s1,s2,n):
        name=n
        nameU=n
        intervalPlot=INTERVAL
        lPlots=[]
        images=[]
        start_time_total=time.time()
        cu=CompU(self.A,self.Er)
        sample=SampPCA(self.A,self.Er)
        ORS=self.Theta
        ORS_PCA=self.Theta
        #SRS=[SampPCA.getISpoint(self.Theta)]
        SRS=SampPCA.getRandIS(self.Theta)
        SRS2=[self.Theta]
        RS=self.Theta
        U=cu.computeUI_Interval(ORS)
        U_PCA=cu.computeU_PCA(ORS_PCA,self.getBasis(SRS))
        t=1
        print()
        print(n)
        print("-----------------\n\n")

        (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
        (X2,Y2)=Visualization(s1,s2,ORS).getPlotsLineFine()
        (X4,Y4)=Visualization(s1,s2,ORS_PCA).getPlotsLineFine()
        lst=SampPCA.getPlotsLineFine(s1,s2,SRS)
        (X3,Y3)=Visualization(s1,s2,SRS2[0]).getPlotsLineFine()
        lPlots=[lst,[(X3,Y3)],[X1,Y1,X2,Y2,X4,Y4]]
        images.append(Visualization.getPlotPCA(s1,s2,lPlots,name+"_0"))

        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress (Optimization): "+str((t*100)/self.T)+"%")
            sys.stdout.flush()

            RS=CompU.prodMatStars(self.A,RS)

            SRS=sample.prodMatIS(SRS)
            SRS2=sample.prodMatStars(SRS2)

            ORS=CompU.addStars(CompU.prodMatStars(self.Ac,ORS),U)
            U=cu.computeUI_Interval(ORS)

            ORS_PCA=CompU.addStars(CompU.prodMatStars(self.Ac,ORS_PCA),U_PCA)
            U_PCA=cu.computeU_PCA(ORS_PCA,self.getBasis(SRS))

            if t%intervalPlot==0:
                (X1,Y1)=Visualization(s1,s2,RS).getPlotsLineFine()
                (X2,Y2)=Visualization(s1,s2,ORS).getPlotsLineFine()
                (X3,Y3)=Visualization(s1,s2,ORS_PCA).getPlotsLineFine()
                lst=SampPCA.getPlotsLineFine(s1,s2,SRS)
                lst2=SampPCA.getPlotsLineFine2(s1,s2,SRS2)
                lPlots=[lst,lst2,[X1,Y1,X2,Y2,X3,Y3]]
                name=n+"_"+str(t)
                images.append(Visualization.getPlotPCA(s1,s2,lPlots,name))

            t=t+1

        print("\n")
        images[0].save("GIFs/"+n+'.gif',save_all=True, append_images=images[1:], optimize=False, duration=250, loop=0)
        print("Total Time: ",time.time()-start_time_total)

    def getBasis(self,samp):
        '''
        Perform PCA of the samples to get a set of
        Basis Vectors.
        Ref: https://www.youtube.com/watch?v=rng04VJxUt4
        '''

        m=len(samp)
        X=np.zeros((m,self.n))

        for i in range(m):
            X[i]=samp[i].reshape((1,self.n))

        Sigma=(1/m)*np.matmul((np.transpose(X)),X)

        (U,S,Vh)=LA.svd(Sigma)

        return U
