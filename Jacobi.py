# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 19:22:28 2018

@author: Youhan0
"""

import numpy as np
from numpy.linalg import * 

def Jacobi(A,b,tol):
    n=len(A)
    D=np.mat(np.zeros((n,n)))
    L=np.mat(np.zeros((n,n)))
    U=np.mat(np.zeros((n,n)))
    for i in range(0,n):
        D[i,i]=A[i,i]
    for i in range(1,n):
        for j in range(0,i):
            L[i,j]=-A[i,j]
    for i in range(0,n):
        for j in range(i+1,n):
            U[i,j]=-A[i,j]
    J=D.I*(L+U)
    G=(D-L).I*U
    w,v=eig(J)
    #print(eig(J))
    #print(eig(G))
    if max(abs(w))>=1:
        print("谱半径大于等于1，不收敛")
        return(0)
    else:
        print("谱半径小于1，收敛")

    N=np.shape(A)[0]
    x=np.mat(np.zeros((N,1)))
    print("初始值：x=")
    print(x)
    def fun(x1):
        x2=np.mat(np.zeros((N,1)))
        for i in range(N):
            sum=0.0
            for j in range(N):
                if j!=i:
                    sum+=A[i,j]*x1[j,0]
            x2[i]=(b[i]-sum)/A[i,i]
        return(x2)
    
    for k in range(50):
        temp=fun(x)
        if max(abs(x[:,0]-temp[:,0]))<tol:
            print("最终结果：x=")
            print(temp) 
            return(temp) 
        print("第" + str(k+1) + "次结果：x=")
        x=temp
        print(x)
    
  
"""main"""
#tol=1e-6 #tol为精度要求
#A=np.mat([[8.0,-3.0,2.0],
#          [4.0,11.0,-1.0],
#          [6.0,3.0,12.0]])
#b=np.mat([[20.0],
#          [33.0],
#          [36.0]])

tol=1e-5 #tol为精度要求
A1=np.mat([[1.0,-5.0,-1.0],
          [4.0,1.0,-1.0],
          [2.0,-1.0,-6.0]])
b1=np.mat([[-8.0],
          [13.0],
          [-2.0]])

tol=1e-5 #tol为精度要求
A2=np.mat([[-2.0,1.0,5.0],
          [4.0,-8.0,1.0],
          [4.0,-1.0,-1.0]])
b2=np.mat([[15.0],
          [-21.0],
          [7.0]])      

Jacobi(A1,b2,tol)
Jacobi(A1,b2,tol)

