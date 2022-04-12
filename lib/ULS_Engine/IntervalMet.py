'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using Interval Arithmetic

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import mpmath as mp
import sys
import time

class Interval:
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
        Indicating array[i][j] has a perturbation of -100*(1-a)% to 100*(b-1)%
        For eg: [0.9,1,1] means a perturbation of -10% and +10%
        '''
        self.Theta=theta # The initial set
        self.T=T #Max number of steps
        self.n=A.shape[0]
        self.A_tilde=self.computeUncertainMat()
        self.methodName="Interval Arithmetic"

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
            A_tilde[key[0]][key[1]]=mp.mpi(min(a,b),max(a,b))
            #A_tilde[key[0]][key[1]]=mp.mpi(self.Er[key][0]*self.A[key[0]][key[1]],self.Er[key][1]*self.A[key[0]][key[1]])

        #print(A_tilde)
        return A_tilde

    def getReachableSetPertFree(self):
        '''
        Computes the reachable set without perturbation
        '''
        At=LA.matrix_power(self.A,self.T)
        rs=np.matmul(At,self.Theta)
        return rs

    def printReachableSetPertFree(self):
        '''
        Computes the reachable set without perturbation
        '''
        start_time=time.time()
        At=LA.matrix_power(self.A,self.T)
        rs=np.matmul(At,self.Theta)
        time_taken=time.time()-start_time
        print()
        print("\n-------------Reachable Set of the Un-perturbed System-------------")
        print("Time Taken: ",time_taken)
        print()
        print(rs)
        print()
        print("Step: ",self.T)
        print("---------------------------------------------------------------")

    def getReachableSet(self):
        '''
        Computes the reachable set of the perturbed system
        '''
        At=LA.matrix_power(self.A_tilde,self.T)
        rs=np.matmul(At,self.Theta)
        return rs

    def printReachableSet(self):
        '''
        Computes the reachable set of the perturbed system
        '''
        start_time=time.time()
        At=LA.matrix_power(self.A_tilde,self.T)
        rs=np.matmul(At,self.Theta)
        time_taken=time.time()-start_time
        print()
        print("\n-------------Reachable Set of the Perturbed System using Interval Arithmetic-------------")
        print("Time Taken: ",time_taken)
        print(rs)
        print()
        print("Step: ",self.T)
        print("---------------------------------------------------------------")



if False:
    A=np.array([
    [1,1,-2],
    [2,0.2,0],
    [0,0.1,0.1]
    ])
    E={
    (0,1):[0.9,1.5],
    (1,0):[0.8,1.2],
    }
    T=2000
    IS=np.array([
    [1],
    [1],
    [1]
    ])
    Interval(A,E,IS,T).printReachableSet()
