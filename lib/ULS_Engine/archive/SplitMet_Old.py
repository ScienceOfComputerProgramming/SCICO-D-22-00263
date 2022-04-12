'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.

Note that this is not the correct version, it uses interval arithmetic.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import mpmath as mp
import sys
import time

class Split:
    '''
    Computes the reachable set of a given linear discrete dynamical system A
    form the initial set Theta, upto time T
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
        self.A_tilde=self.computeUncertainMat()
        self.Ac=self.computeCenter()
        self.methodName="Split"

    def computeCenter(self):
        '''
        Computes the center point matrix of the interval matrix self.A
        '''
        Ac=np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                if (isinstance(self.A_tilde[i][j],int)) or (isinstance(self.A_tilde[i][j],float)) or (isinstance(self.A_tilde[i][j],np.float128)):
                    Ac[i][j]=self.A_tilde[i][j]
                else:
                    '''print(mp.nstr(self.A_tilde[i][j]))
                    print(mp.nstr(self.A_tilde[i][j]).split(',')[0][1:])
                    print(mp.nstr(self.A_tilde[i][j]).split(',')[1][:-1])'''
                    a=float(mp.nstr(self.A_tilde[i][j]).split(',')[0][1:])
                    b=float(mp.nstr(self.A_tilde[i][j]).split(',')[1][:-1])
                    c=(a+b)/2
                    Ac[i][j]=c
        #print(Ac)
        return Ac

    def computeUncertainMat(self):
        '''
        Computes the interval uncertain matrix
        '''
        A_tilde=np.zeros((self.n,self.n),dtype=object)
        for i in range(self.n):
            for j in range(self.n):
                A_tilde[i][j]=self.A[i][j]
        for key in self.Er:
            a=float(self.Er[key][0]*self.A[key[0]][key[1]])
            b=float(self.Er[key][1]*self.A[key[0]][key[1]])
            #print(a,b)
            #print(mp.mpi(float(a),float(b)))
            #A_tilde[key[0]][key[1]]=mp.mpi(float(self.Er[key][0]*self.A[key[0]][key[1]]),float(self.Er[key][1]*self.A[key[0]][key[1]]))
            #A_tilde[key[0]][key[1]]=mp.mpi(self.Er[key][0]*self.A[key[0]][key[1]],self.Er[key][1]*self.A[key[0]][key[1]])
            A_tilde[key[0]][key[1]]=mp.mpi(min(a,b),max(a,b))

        #print(A_tilde)
        return A_tilde


    def computeU(self,rs):
        '''
        Computes the effect of uncertainty on the reachable set
        '''
        #A_tildeX=np.matmul(self.A_tilde,rs)
        #AcX=np.matmul(self.Ac,rs)
        #U=A_tildeX-AcX
        diff=self.A_tilde-self.Ac
        U=np.matmul(diff,rs)
        return U

    def getReachableSet(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        ORS=self.Theta
        U=self.computeU(ORS)
        t=1
        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Progress: "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            ORS=np.matmul(self.Ac,ORS)+U
            U=self.computeU(ORS)
            t=t+1
        print()
        return ORS

    def printReachableSet(self):
        '''
        Implements the main algorithm of splitting the effect of the constant and
        the uncertain part
        '''
        start_time=time.time()
        ORS=self.Theta
        U=self.computeU(ORS)
        t=1
        print()
        while (t<=self.T):
            sys.stdout.write('\r')
            sys.stdout.write("Splitting Algorithm Progress: "+str((t*100)/self.T)+"%")
            sys.stdout.flush()
            ORS=np.matmul(self.Ac,ORS)+U
            U=self.computeU(ORS)
            t=t+1
        time_taken=time.time()-start_time
        print()
        print("\n-------------Reachable Set of the Perturbed System using Splitting Method-------------")
        print("Time Taken: ",time_taken)
        print(ORS)
        print("---------------------------------------------------------------")

if False:
    A2=np.array([
    [1,1,-2],
    [2,0.2,0],
    [0,0.1,0.1]
    ])
    A=np.array([
    [2,1,-2,1],
    [2,1.2,0,0],
    [0,0.1,1.1,1],
    [0,0,0,1]
    ])
    E={
    (0,1):[0.9,1.5],
    (1,0):[0.8,1.2],
    (2,3):[0.9,1.1]
    }
    E2={
    (0,1):[0.9,1.5],
    (1,0):[0.8,1.2],
    }
    T=2000
    IS2=np.array([
    [1],
    [1],
    [1]
    ])
    IS=np.array([
    [1],
    [1],
    [1],
    [1]
    ])
    Split(A,E,IS,T).printReachableSet()
