'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the aprroximate
reachable set of the uncertain system using random sampling of the
uncertain system

Documentation: Not yet available. (TODO)
'''

import random
import math
import numpy as np
from ComputeU import *

from VisualizationReachSet import *

PRECISION=1e25

NO_SAMPPLES=50


class SampPCA:
    '''
    Given an uncertain linear system, this class provides the
    necessary functions to sample the system at random points
    and compute its reachable set
    '''

    def __init__(self,A,er,i=NO_SAMPPLES):
        self.A=A
        self.Er=er
        self.n=self.A.shape[0]
        self.nSamples=i
        self.samples=self.getSamplesULS()

    def getSamplesULS(self):
        '''
        Returns random samples of the uncertain linear
        system.
        '''
        sam=[]

        for i in range(self.nSamples):
            sam.append(self.getRandMat())

        '''for s in sam:
            print(s[3][7],s[4][6])
            print()
        exit(0)'''

        return sam

    def getRandIS(starSet):
        is_list=[]
        for i in range(NO_SAMPPLES):
            is_list.append(SampPCA.getISpoint(starSet))
        return is_list

    @staticmethod
    def getISpoint(starSet):
        C=starSet[0]
        V=starSet[1]
        P=starSet[2]
        (n,m)=V.shape

        P_point=[]
        IS=np.zeros((n,1))
        for p in P:
            lb=p[0]*PRECISION
            ub=p[1]*PRECISION
            if lb!=ub:
                rd=random.randrange(lb,ub)
                rd=rd/PRECISION
            else:
                rd=lb/PRECISION
            P_point.append(rd)

        for i in range(n):
            obj=0
            for j in range(m):
                obj=obj+(P_point[j]*V[i][j])
            IS[i][0]=C[i]+obj

        #print(IS)

        return IS

    @staticmethod
    def getISpoint2(starSet):
        C=starSet[0]
        V=starSet[1]
        P=starSet[2]
        (n,m)=V.shape

        P_point=[]
        IS=np.zeros((n,1))
        for p in P:
            lb=p[0]
            ub=p[1]
            if lb!=ub:
                rd=(lb+ub)/2
                rd=rd
            else:
                rd=lb
            P_point.append(rd)

        for i in range(n):
            obj=0
            for j in range(m):
                obj=obj+(P_point[j]*V[i][j])
            IS[i][0]=C[i]+obj

        return IS

    def getRandMat(self):
        '''
        Returns a random matrix of the uncertain linear
        system.
        '''

        rMat=np.zeros((self.n,self.n),dtype=float)

        for i in range(self.n):
            for j in range(self.n):
                if (i,j) in self.Er:
                    a=self.A[i][j]*self.Er[(i,j)][0]
                    b=self.A[i][j]*self.Er[(i,j)][1]
                    lb=min(a,b)
                    ub=max(a,b)
                    #print(self.A[i][j])
                    #print(lb,ub)
                    #print(math.ceil(lb*PRECISION),math.floor(ub*PRECISION))
                    if lb==0 and ub==0:
                        r=0
                    else:
                        rd=random.randrange(math.ceil(lb*PRECISION),math.floor(ub*PRECISION))
                        r=float(rd/PRECISION)
                        if (r>ub or r<lb):
                            print("Gotcha!!")
                            print("r :",r)
                            print("[lb,ub]: ",lb,ub)
                            exit(0)
                    rMat[i][j]=r
                else:
                    rMat[i][j]=self.A[i][j]

        #print(rMat[3][7],rMat[4][6])
        #exit(0)

        return rMat

    def prodMatIS(self,RS_list):
        '''
        Multiply the random samples with RS_list
        '''

        is_list=[]

        if (len(RS_list)==1):
            for s in self.samples:
                is_list.append(np.matmul(s,RS_list[0]))
        else:
            for i in range(len(RS_list)):
                is_list.append(np.matmul(self.samples[i],RS_list[i]))


        '''if (len(RS_list)==1):
            for s in self.samples:
                is_list.append(np.matmul(self.A,RS_list[0]))
        else:
            for i in range(len(RS_list)):
                is_list.append(np.matmul(self.A,RS_list[i]))'''

        return is_list

    def prodMatStars(self,RS_list):
        '''
        Multiply the random samples with RS_list
        '''
        star_list=[]

        if (len(RS_list)==1):
            for s in self.samples:
                star_list.append(CompU.prodMatStars(s,RS_list[0]))
        else:
            for i in range(len(RS_list)):
                star_list.append(CompU.prodMatStars(self.samples[i],RS_list[i]))

        return star_list

    @staticmethod
    def getPlotsLineFine(th1,th2,RS_list):
        '''
        Returns the plots for the given RS_list
        '''

        plots=[]

        for s in RS_list:
            X=s[th1]
            Y=s[th2]
            plots.append((X,Y))

        return plots

    @staticmethod
    def getPlotsLineFine2(th1,th2,RS_list):
        '''
        Returns the plots for the given RS_list
        '''

        plots=[]

        for star in RS_list:
            (X,Y)=Visualization(th1,th2,star).getPlotsLineFine()
            plots.append((X,Y))

        return plots



if False:
    A=np.array([
    [2,1],
    [2,1.2]
    ])
    E={
    (0,1): [0.98,1.02]
    }
    s=Sampling(A,E)
    C=[0,0]
    V=np.array([
    [1,0],
    [0,1]
    ])
    P=[(-1,1),(-1,1)]
    rs=(C,V,P)
    print(Sampling.getPlotsLineFine(0,1,s.prodMatStars([rs])))
