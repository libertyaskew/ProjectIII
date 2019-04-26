"""
RATIONAL FUNCTION RECONSTRUCTION

Thiele implements the rational function reconstruction algorithm for the black box function f 
as presented in the paper.

Code for plots in paper is in Plot_Time and Plot_MaxInt and Plot_Accuracy
"""

import numpy as np
from random import * 
from sympy import *
from fractions import Fraction
import time
from matplotlib import pyplot as plt

def g(x):
    return (Fraction( (3*x**2+2*x+ 2),(x-6)))

def Thiele(g , ran):
    #returns Thiele coefficents
    #Thiele that evaluates function until  reciprocal difference is 0
        #if repeated integers sampled, algorithm fails
        x = []
        while (len(x)) < 4:
            while True:
                R = randint(1,ran)
                try:
                    g(R)
                    if R not in x:
                        x = np.append(x,R)
                        break
                except ZeroDivisionError:
                    continue
        x = [int(x[i]) for i in range (4)]
        yn = np.array([g(x[i]) for i in range (4)])
        xn = x.copy()
        c = []
        while True:
            m = len(xn)
            a = yn.copy()
            b = []
            for k in range (1,m):
                b.append(a[k:m-1].copy())
                if a[-1]==a[-2]:
                    return(c[-1], xn)
                a[k:m] = [ Fraction((xn[i+k]-xn[i]),(a[i+k] - a[i+k-1])) for i in range (m-k)]
                if k>1:
                    a[k:m]+=list(b[k-2])
            for k in range (2 , len(a)):
                a[k]-= a[k-2]
            while True:
                R = randint(1,ran)
                try:
                    g(R)
                    if R not in xn:
                        xn = np.append(xn,R)
                        break
                except ZeroDivisionError:
                    continue
                else:
                    continue
            yn = np.append(yn , g(xn[-1]))
            c.append(a)
            
def Plot2(): 
    #produces plot for each #S of time taken
    S = [20+5*i for i in range (17)]  
    tl = [] #time list
    for j in range (len(S)):
        c = i = tc = 0 #tc = time counter, c = error counter , i = iteration counter
        for i in range (10**3):
            start = time.time()
            p = Thiele3(g, S[j])
            elapsed = time.time()-start
            tc += elapsed
            i += 1
        tl.append(tc/10**3)
        print(tl)
    plt.plot(S,tl,'x')
    plt.show()
Plot2()

def Test(f,n):
    #reconstructs rational function endless number of times to test accuracy 
    i=0
    while True:
        T = (Thiele3(f,50))
        i +=1
        if i%10000==0:
            print (i) #visual counter to show it is working
        if len(T[1]) != n: #adjust to be expected length of xn for correct reconstruction
            print (T, i , 'ERROR') #visual display of error
