'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com
- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.
Documentation: Not yet available. (TODO)
'''

import math
import numpy as np
import copy
import mpmath as mp

class Grid:
    '''
    Grids the given inital set into pieces
    '''

    def __init__(self,Theta,st1,st2):

        self.Theta=Theta
        self.state1=st1
        self.state2=st2

    def splitStar(self):
        '''
        Performs gridding just by dividing the predicate.
        -6<=alpha_i<=6 => -6<=alpha_i<=0, 0<=alpha_i<=-6
        '''
        C=self.Theta[0]
        V=self.Theta[1]
        P=self.Theta[2]

        l1=abs(P[self.state1][0]-P[self.state1][1])
        mid1=math.ceil(l1/2)
        l2=abs(P[self.state2][0]-P[self.state2][1])
        mid2=math.ceil(l2/2)

        p1=(P[self.state1][0],P[self.state1][0]+mid1)
        p2=(P[self.state2][0],P[self.state2][0]+mid2)
        p3=(P[self.state1][0]+mid1,P[self.state1][1])
        p4=(P[self.state2][0]+mid2,P[self.state2][1])

        P1=copy.copy(P)
        P2=copy.copy(P)
        P3=copy.copy(P)
        P4=copy.copy(P)

        P1[self.state1]=p1
        P1[self.state2]=p2
        P2[self.state1]=p1
        P2[self.state2]=p4
        P3[self.state1]=p3
        P3[self.state2]=p2
        P4[self.state1]=p3
        P4[self.state2]=p4

        listStars=[(C,V,P1),(C,V,P2),(C,V,P3),(C,V,P4)]

        return listStars

    def printStars(listStars):

        for star in listStars:
            print("---Star---")
            print("Center: ",star[0])
            print("Vector: \n",star[1])
            print("Predicate: ",star[2])

class GridMat:
    '''
    Grids the given uncertain matrix into pieces
    '''

    def __init__(self,A,Er):

        self.A=A
        self.Er=Er

    def splitEr(self):
        '''
        Breaks the interval matrix to 4 interval matrices
        by spiltting two random intervals from state1
        and state2
        '''

        n=self.A.shape[0]
        Er1=copy.copy(self.Er)
        Er2=copy.copy(self.Er)
        Er3=copy.copy(self.Er)
        Er4=copy.copy(self.Er)

        l=len(self.Er)

        if l==1:
            (x1,y1)=list(self.Er.keys())[0]
            mid1=(self.Er[(x1,y1)][1]-self.Er[(x1,y1)][0])/2
            I1=(self.Er[(x1,y1)][0],self.Er[(x1,y1)][0]+mid1)
            I2=(self.Er[(x1,y1)][0]+mid1,self.Er[(x1,y1)][1])
            Er1[(x1,y1)]=I1
            Er2[(x1,y1)]=I2

            return [Er1,Er2,Er1,Er2]

        else:
            (x1,y1)=list(self.Er.keys())[0]
            (x2,y2)=list(self.Er.keys())[l-1]
            mid1=(self.Er[(x1,y1)][1]-self.Er[(x1,y1)][0])/2
            mid2=(self.Er[(x2,y2)][1]-self.Er[(x2,y2)][0])/2
            I1=[self.Er[(x1,y1)][0],self.Er[(x1,y1)][0]+mid1]
            I2=[self.Er[(x1,y1)][0]+mid1,self.Er[(x1,y1)][1]]
            I3=[self.Er[(x2,y2)][0],self.Er[(x2,y2)][0]+mid2]
            I4=[self.Er[(x2,y2)][0]+mid2,self.Er[(x2,y2)][1]]
            Er1[(x1,y1)]=I1
            Er1[(x2,y2)]=I3
            Er2[(x1,y1)]=I1
            Er2[(x2,y2)]=I4
            Er3[(x1,y1)]=I2
            Er3[(x2,y2)]=I3
            Er4[(x1,y1)]=I2
            Er4[(x2,y2)]=I4

            erList=[Er1,Er2,Er3,Er4]

            return erList

    def splitErOne(self):
        '''
        Breaks the interval matrix to 2 interval matrices
        by spiltting two random intervals from state1
        and state2
        '''

        n=self.A.shape[0]
        Er1={}
        Er2={}

        l=len(self.Er)

        (x1,y1)=list(self.Er.keys())[0]
        mid1=(self.Er[(x1,y1)][1]-self.Er[(x1,y1)][0])/2
        I1=(self.Er[(x1,y1)][0],self.Er[(x1,y1)][0]+mid1)
        I2=(self.Er[(x1,y1)][0]+mid1,self.Er[(x1,y1)][1])
        Er1[(x1,y1)]=I1
        Er2[(x1,y1)]=I2

        return [Er1,Er2]

    def splitErOne2(self):
        '''
        Breaks the interval matrix to 2 interval matrices
        by spiltting two random intervals from state1
        and state2
        '''

        n=self.A.shape[0]
        Er1={}
        Er2={}
        Er3={}
        Er4={}


        (x1,y1)=list(self.Er.keys())[0]
        samp=(self.Er[(x1,y1)][1]-self.Er[(x1,y1)][0])/4
        I1=(self.Er[(x1,y1)][0],self.Er[(x1,y1)][0]+samp)
        I2=(self.Er[(x1,y1)][0]+samp,self.Er[(x1,y1)][0]+(2*samp))
        I3=(self.Er[(x1,y1)][0]+(2*samp),self.Er[(x1,y1)][0]+(3*samp))
        I4=(self.Er[(x1,y1)][0]+(3*samp),self.Er[(x1,y1)][1])
        Er1[(x1,y1)]=I1
        Er2[(x1,y1)]=I2
        Er3[(x1,y1)]=I3
        Er4[(x1,y1)]=I4

        return [Er1,Er2,Er3,Er4]


    def printEr(listMat):

        for mat in listMat:
            print("--Error--")
            print(mat)

if False:
    C=[0,0,0]
    V=np.array([
    [1,-1,1],
    [-1,1,0],
    [1,0,-1],
    ])
    P=[(-6,5),(-7,9),(-2,9)]
    rs=(C,V,P)
    Grid.printStars(Grid(rs,0,2).splitStar())

if False:
    A=np.array([
    [2,1,-2,1],
    [2,1.2,0,0],
    [0,0.1,1.1,1],
    [0,0,0,1]
    ])
    E={
    (0,1):[0.9,1.1],
    (1,0):[0.8,1.2],
    (2,3):[0.9,1.1]
    }
    GridMat.printEr(GridMat(A,E).splitEr())
