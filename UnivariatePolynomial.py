"""
UNIVARIATE POLYNOMIAL

Codes for univariate black box reconstruction and plots for testing accuracy and time.
"""
import numpy as np
from sympy import *
from random import *
from matplotlib import pyplot as plt
import time

def PolyNewton_2(f , ran):
#black box reconstruction for a function f and a fixed size of maximum S
        x = np.array([randint(1,ran)])
        while (len(x)) < 3:
            R = randint(1,ran)
            if R not in x:
                x = np.append(x,R)
        yn = np.array([f(x[i]) for i in range (3)])
        xn = x.copy()
        while True:
            m = len (xn)
            a = yn.copy()
            for k in range (1,m):
                a[k:m] = (a[k:m] - a[k-1:m-1])/(xn[k:m] - xn[0:m-k])
            while True:
                R = randint(1,ran)
                if R not in xn:
                    xn = np.append(xn,R)
                    break
                else:
                    continue
            yn = np.append(yn , f(xn[-1]))
            if a[-1]==0:
                break
        xl = list(xn)
        poly2 = 0
        for i in range( len (a)):
            d = int(a[i])
            poly2 += d * np.poly1d(xl[:i],True)
        #returns number of iterations, reconstructed polynomial and the list of xn values sampled
        #this information is needed for later functions
        return (len(a) , poly2, xn)

    
def Plot_Failure(f):
    #produces plot for each #S probability of failure 
    S = [10+5*i for i in range (19)]
    l = []  #list of probabilities for each S
    for j in range (len (S)):
        i = 0
        c = 0 #tc = time counter, c = error counter , i = iteration counter
        N = 10**3
        for i in range (N):
            p = PolyNewton_2(f, S[j])
            i += 1
            if p[0]!= 6:
                c += 1
        l.append(c/N)
    plt.plot(S,l,'x')
    plt.xlabel('#S')
    plt.ylabel('Probability of incorrect reconstruction')
    plt.show()

  
def Plot_Time(f): 
    #produces plot for each #S of time taken
    S = [20+100*i for i in range (30)]  
    tl = [] #time list
    for j in range (len(S)):
        c = i = tc = 0 #tc = time counter, c = error counter , i = iteration counter
        for i in range (10**3):
            start = time.time()
            p = PolyNewton_2(f, S[j])
            elapsed = time.time()-start
            tc += elapsed
            i += 1
        tl.append(tc/10**3)
    plt.plot(S,tl,'x')
    plt.xlabel('#S')
    plt.ylabel('Time taken for reconstruction')
    plt.show()
   
    
def Test(f , n): 
    i = c = 0
    # finds probability of failure for fixed S for endless iterations for function f
    #set n to be the number of iterations required for correct reconstruction
    while True:
        p = PolyNewton_2(f, 100)
        i += 1 #iteration counter
        if i%1e7 == 0:
            print (i) #visual counter
        if p[0]!= n:  # number of reconstructions indicates early termination and failure
            print (p , c , i) #visual display of error
            c += 1 #error counter


