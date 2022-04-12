'''
Author: Bineet Ghosh, under supervision of Dr. Sridhar Duggirala
Email: ghosh.bineet22@gmail.com

- Given a Linear System with uncertainties, compute the reachable set of the
uncertain system using the Splitting Method.

Documentation: Not yet available. (TODO)
'''

import numpy as np
import numpy.linalg as LA
import mpmath as mp
import math
import time

PRINT_COUNT=5
EPSILON=1e-2
P=1

class OrdUnc:
    '''
    Orders the cells of the matrix.
    Based on the paper:
    "Perturbation Results for Singular Values"
    by Torsten Soderstrom.

    Cells which affects the singular values
    more than others are more sensitive to
    perturbation.
    '''

    def __init__(self,A):

        self.A=A
        self.n=A.shape[0]

    def getOrder(self):

        '''
        Returns an ordering of the cells based on
        decreasing sensitivity to perturbation
        '''
        ordMat=self.getSVSentivity()
        ord=OrdUnc.sortMat(ordMat)
        return ord

    def getSVSentivity(self,pert=1):
        '''
        Run the ordering algorithm
        '''
        ordMat=np.zeros((self.n,self.n))
        B=np.zeros((self.n,self.n))

        if self.determineCase()==1:
            for i in range(self.n):
                for j in range(self.n):
                    B[i][j]=pert
                    # Processing
                    ordMat[i][j]=self.distinctPos(B)
                    # End: Processing
                    B[i][j]=0
        else:
            for i in range(self.n):
                for j in range(self.n):
                    #print("(",i,",",j,")")
                    B[i][j]=pert
                    # Processing
                    ordMat[i][j]=self.multSig(B)
                    # End: Processing
                    B[i][j]=0
                    #print("--\n")
        return ordMat

    def getOrderRelative(self):

        '''
        Returns an ordering of the cells based on
        decreasing sensitivity to perturbation
        '''

        ordMat=np.zeros((self.n,self.n))
        B=np.zeros((self.n,self.n))

        if self.determineCase()==1:
            for i in range(self.n):
                for j in range(self.n):
                    B[i][j]=self.A[i][j]
                    # Processing
                    ordMat[i][j]=self.distinctPos(B)
                    # End: Processing
                    B[i][j]=0
        else:
            for i in range(self.n):
                for j in range(self.n):
                    B[i][j]=self.A[i][j]
                    # Processing
                    ordMat[i][j]=self.multSig(B)
                    # End: Processing
                    B[i][j]=0

        #print(ordMat)
        ord=OrdUnc.sortMat(ordMat)

        return ord

    def determineCase(self):
        '''
        Determine which equation to follow
        '''

        distinct=False
        zero=False
        allPositive=True
        (u,sVals,vh)=LA.svd(self.A)

        for s in sVals:
            if s<=0:
                zero=True
                allPositive=False
                break
        if len(sVals)==len(set(sVals)):
            distinct=True

        if distinct==True and allPositive==True:
            return 1
        elif distinct==False and allPositive==True:
            return 2
        else:
            return 3

    def distinctPos(self,B):

        '''
        Case 1: Distinct and All Positive
        '''

        (u,sVals,vh)=LA.svd(self.A)
        uh=u.conjugate().transpose()
        v=vh.conjugate().transpose()
        #print(v)
        #print(v[0,:])
        #print(uh[:,0].reshape(self.n,1))
        dSig=[]

        max=-9999
        prodQ=1
        for j in range(self.n):
            uj=uh[:,j]
            vj=v[j,:].reshape(self.n,1)
            q=np.matmul(np.matmul(uj,B),vj)
            prodQ=prodQ*q
            if q.real[0]>max:
                max=q.real[0]

        return prodQ

    def multSig(self,B):

        '''
        Case 2 and 3: Multiple
        '''

        (u,sVals,vh)=LA.svd(self.A)
        '''print(u)
        print(sVals)
        print(vh)'''
        l={}

        start=0
        end=0
        for i in range(self.n):
            if sVals[i]!=sVals[start]:
                end=i-1
                sig=sVals[start]
                uSig=(u[:,start:(end+1)])
                vSig=(vh[start:(end+1),:])
                l[sig]=[uSig,vSig]
                start=i
                end=i

        sig=sVals[start]
        uSig=(u[:,start:(end+2)])
        vSig=(vh[start:(end+2),:])
        l[sig]=[uSig,vSig]

        '''print("\n-")
        for keys in l:
            print("Sigma: ",keys)
            print(l[keys][0])
            print(l[keys][1])
            print("-")
        '''

        max=-9999
        for keys in l:
            sig=keys
            uSig=l[keys][0]
            vSig=l[keys][0]
            uSigH=uSig.conjugate().transpose()
            vSigH=vSig.conjugate().transpose()
            Bh=B.conjugate().transpose()
            prodQ=1
            if sig>0:
                T=(np.matmul(np.matmul(vSigH,Bh),uSig))+(np.matmul(np.matmul(uSigH,B),vSig))
                (eVals,w)=LA.eig(T)
                for eig in eVals:
                    q=sig+((EPSILON*(eig))/2)+(EPSILON*EPSILON)
                    prodQ=prodQ*q
                    #print("\tq: ",q)
                    if q>max:
                        max=q
            else:
                #vSigH*Bh*uSig*uSigH*B*vSig
                T=np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(vSigH,Bh),uSig),uSigH),B),vSig)
                (eVals,w)=LA.eig(T)
                for eig in eVals:
                    q=EPSILON*(math.sqrt(eig))+(EPSILON*EPSILON)
                    prodQ=prodQ*q
                    #print("\tq: ",q)
                    if q>max:
                        max=q

        return prodQ

    def sortMat(A):

        '''
        Sorts the matrix
        '''
        n=A.shape[0]
        lst=[]
        for i in range(n):
            for j in range(n):
                max=A[i][j]
                indx=(i,j)
                for k in range(n):
                    for l in range(n):
                        if A[k][l]>max:
                            indx=(k,l)
                            max=A[k][l]
                            #print(max)
                lst.append(indx)
                A[indx[0]][indx[1]]=-9999

        return lst

    def printReport(self):

        '''
        Prints the cells in decreasing order of
        sensitivity towards perturbation
        '''

        start_time=time.time()
        ord=self.getOrder()
        end_time=time.time()-start_time
        #print(ord)

        ct=1
        print("----- Sensitivity of the cells in decreasing order -----")
        for o in ord:
            print(o)
            ct=ct+1
            if ct>PRINT_COUNT:
                print("...",(self.n*self.n-PRINT_COUNT)," more ...")
                break
        print("Time Taken: ",end_time)
        c=self.determineCase()
        if (c==1):
            print("Case: Distinct and Positive")
        elif (c==2):
            print("Case: Multiple and Positive")
            (u,sVals,vh)=LA.svd(self.A)
            print("Singular Values: ",sVals)
        else:
            print("Case: Zero")
            (u,sVals,vh)=LA.svd(self.A)
            print("Singular Values: ",sVals)
        print("---------------------<END>--------------------")

    def printReportCompare(self):

        '''
        Prints the cells in decreasing order of
        sensitivity towards perturbation
        '''

        start_time=time.time()
        ord=self.getOrder()
        ordRel=self.getOrderRelative()
        end_time=time.time()-start_time
        l=len(ord)
        #print(ord)

        ct=1
        print("Relative == Absolute: ",ord==ordRel)
        print("Relative == Absolute (top 5): ",ord[:PRINT_COUNT]==ordRel[:PRINT_COUNT])
        print("Relative == Absolute (bottom 5): ",ord[l-PRINT_COUNT:]==ordRel[l-PRINT_COUNT:])
        print("Time Taken: ",end_time)
        print("Details>>>>")
        print("----- Sensitivity of the top cells in decreasing order (Relative)-----")
        for o in ordRel:
            print(o)
            ct=ct+1
            if ct>PRINT_COUNT:
                print("...",(self.n*self.n-PRINT_COUNT)," more ...")
                break

        print("----- Sensitivity of the bottom cells in decreasing order (Relative)-----")
        for o in ordRel[l-PRINT_COUNT:]:
            print(o)
            ct=ct+1
        print("...",(self.n*self.n-PRINT_COUNT)," more ...")

        Er=OrdUnc.dictionarify(ordRel[:PRINT_COUNT])
        Er2=OrdUnc.dictionarify(ordRel[l-PRINT_COUNT:])

        #print(Er)
        #print(Er2)

        return (Er,Er2)

        print("---------------------<END>--------------------")

    def dictionarify(ord):
        Er={}

        for o in ord:
            Er[o]=[1-(P/100),1+(P/100)]

        return Er


if False:
    A=np.array([
    [1,0,0,2],
    [0,5,0,0],
    [3,0,1,0],
    [2,0,0,1],
    ])
    OrdUnc(A).testFun()
