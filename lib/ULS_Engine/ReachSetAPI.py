'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the provided list of methods

Documentation: Not yet available. (TODO)
'''

h=0.01 # Discretization Parameter

from SplitMet import *
from IntervalMet import *

class ReachSet:
    '''
    Computes the reachable set of the given linear uncertain system
    '''

    def __init__(self,A,B,mode,E,theta,T,mList=[]):
        self.dynA=A # Linear system matrix A (without perturbation)
        self.dynB=B # Linear system matrix B (without perturbation)
        self.mode=mode
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
        self.methodList=mList # List of method names based on which the reach set will be computed
        self.A=ReachSet.createMatrix(self.dynA,self.dynB,self.mode,h)

    def getReachableSetPertFree(self):
        '''
        Returns the reachable set of the system without perturbation
        '''
        return Interval(self.A,self.Er,self.Theta,self.T).getReachableSetPertFree()

    def printReachableSetPertFree(self):
        '''
        Prints the reachable set of the system without perturbation
        '''
        Interval(self.A,self.Er,self.Theta,self.T).printReachableSetPertFree()

    def getReachableSet(self,method):
        '''
        Returns the reachable set of the perturbed system based on the given method
        '''
        if method.lower()=='split':
            return Split(self.A,self.Er,self.Theta,self.T).getReachableSet()
        elif method.lower()=='interval':
            return Interval(self.A,self.Er,self.Theta,self.T).getReachableSet()

    def printReachableSet(self,method):
        '''
        Returns the reachable set of the perturbed system based on the given method
        '''
        if method.lower()=='split':
            Split(self.A,self.Er,self.Theta,self.T).printReachableSet()
        elif method.lower()=='interval':
            Interval(self.A,self.Er,self.Theta,self.T).printReachableSet()

    def getReachableSetAll(self):
        '''
        Returns a list of reachable sets computed based
        on the methods provided in the list
        '''
        rsList=[]
        for m in self.methodList:
            if m.lower()=='split':
                rsList.append(Split(self.A,self.Er,self.Theta,self.T).getReachableSet(),m)
            elif m.lower()=='interval':
                rsList.append(Interval(self.A,self.Er,self.Theta,self.T).getReachableSet(),m)

        return rsList

    def printReachableSetAll(self):
        '''
        Prints a list of reachable sets computed based
        on the methods provided in the list
        '''
        for m in self.methodList:
            if m.lower()=='split':
                Split(self.A,self.Er,self.Theta,self.T).printReachableSet()
            elif m.lower()=='interval':
                Interval(self.A,self.Er,self.Theta,self.T).printReachableSet()

    def compareReachSets(self):
        '''
        Compares the reachable set computed based on various approaches
        '''
        for m in self.methodList:
            if m.lower()=='split':
                Split(self.A,self.Er,self.Theta,self.T).printReachableSet()
            elif m.lower()=='interval':
                Interval(self.A,self.Er,self.Theta,self.T).printReachableSet()
        Interval(self.A,self.Er,self.Theta,self.T).printReachableSetPertFree()



    @staticmethod
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
            I=np.zeros((n1,n1),dtype=np.float)
            for i in range(n1):
                I[i][i]=1
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





if False:
    A=np.array([
    [1,1,-2],
    [2,0.2,0],
    [0,0.1,0.1]
    ])
    B=np.array([
    [1],
    [0],
    [1]
    ])
    E={
    (0,1):[0.9,1.1],
    (1,0):[0.8,1.2],
    (2,3):[0.9,1.1]
    }
    mode='.'
    T=2000
    IS=np.array([
    [1],
    [1],
    [1],
    [1]
    ])
    mList=['Split','Interval']
    rs=ReachSet(A,B,mode,E,IS,T,mList)
    rs.compareReachSets()
