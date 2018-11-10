# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 17:27:50 2018

@author: 闫汝海
"""

import sympy as sp

def dichotomy(y,x,a,b,tol):
    
    if (y.subs(x,a)<0 and y.subs(x,b)<0) or (y.subs(x,a)>0 and y.subs(x,b)>0):
        print("区间端点同号！")
        return(0)
    if  y.subs(x,a)==0:
        print("y(a)=0")
        return(0)
    if  y.subs(x,b)==0:
        print("y(b)=0")
        return(0)
    
    if(abs(a-b)>tol/2):
        if  y.subs(x,a)<0 and  y.subs(x,b)>0:
            if y.subs(x,(a+b)/2)==0:
                print("y(("+str(a)+"+"+str(b)+")/2)=0")
                return(y.subs(x,(a+b)/2))
            if y.subs(x,(a+b)/2)<0:
                print("a=" + str((a+b)/2) + ",b=" + str(b))
                dichotomy(y,x,(a+b)/2,b,tol)
            if y.subs(x,(a+b)/2)>0:
                print("a=" + str(a) + ",b=" + str((a+b)/2))
                dichotomy(y,x,a,(a+b)/2,tol)
            
        if y.subs(x,a)>0 and  y.subs(x,b)<0:
            if y.subs(x,(a+b)/2)==0:
                print("y(("+str(a)+"+"+str(b)+")/2)=0")
                return(y.subs(x,(a+b)/2))
            if y.subs(x,(a+b)/2)<0:
                print("a=" + str(a) + ",b=" + str((a+b)/2))
                dichotomy(y,x,a,(a+b)/2,tol)
            if y.subs(x,(a+b)/2)>0:
                print("a=" + str((a+b)/2) + ",b=" + str(b))
                dichotomy(y,x,(a+b)/2,b,tol)
    else:
    #if(abs(a-b)<=tol/2):
        print("近似值为：x=" + str((a+b)/2) + ",y=" + str(y.subs(x,(a+b)/2)))
        return(y.subs(x,(a+b)/2))

    
"""main"""
x=sp.symbols("x")
y=sp.symbols("y")

"""以下是函数以及精度要求"""
y=x**3-x-1
a=1.0
b=1.5
tol=1e-2
#a,b是区间左右端点
#精度为1e-2

dichotomy(y,x,a,b,tol)

