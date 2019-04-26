"""
3 VARIATE POLYNOMIAL RECONSTUCTION

Works by calling Newton_3V(f) where f is a function of 3 variables. Implement 
reconstruction algorithm as presented in paper

Plot_Failure and Plot_Time is the code used to produce the plots 
seen in the paper

"""
import numpy as np
from sympy import *
from fractions import Fraction
import time
from matplotlib import pyplot as plt
from random import *

def PolyNewtonZ (f , xs , yt,M):
#returns a_st(z) for each xs in xn and yt in yn
        zn = np.array([randint(1,M)])
        while (len(zn)) < 3:
            R = randint(1,M)
            if R not in zn:
                zn = np.append(zn,R)
        z  = symbols ('z')
        fn =[f(xs,yt,zn[i]) for i in range (len(zn))]
        while True:
            m = len (zn)
            a = fn.copy()
            for k in range (1,m):
                a[k:m] = [(a[i+k] - a[i+k-1])/(zn[i+k]-zn[i]) for i in range (m-k)]
            while True:
                R = randint(1,M)
                if R not in zn:
                    zn = np.append(zn,R)
                    break
                else:
                    continue
            fn = np.append(fn , f(xs,yt,zn[-1]))
            if a[-1] ==0:
                break
        xl = list(zn)
        poly = 0
        XProd = 1
        for i in range( len (a)):
            d = int(a[i])
            poly += d * XProd
            XProd = XProd * (z - xl[i])
        return (simplify(poly))
    
def PolyNewtonY (yn , zn):
#returns a_s(y,z) using zn = a_st(z) 
        yn = np.array(yn)
        zn = np.array(zn)
        y = symbols('y')
        m = len(zn)
        a = zn.copy()
        for k in range (1,m):
            a[k:m] = [(a[i+k] - a[i+k-1])/(yn[i+k]-yn[i]) for i in range (m-k)]
        xl = list(yn)
        poly = 0
        XProd = 1
        for i in range( len (a)):
            d = a[i]
            poly += d * XProd
            XProd = XProd * (y - xl[i])
        return (simplify(poly))
    
def PolyNewtonX (xn , yn):
#returns f(x,y,z) using yn = a_s(y,z)
        xn = np.array(xn)
        yn = np.array(yn)
        x = symbols('x')
        m = len(yn)
        a = yn.copy()
        for k in range (1,m):
            a[k:m] = [(a[i+k] - a[i+k-1])/(xn[i+k]-xn[i]) for i in range (m-k)]
        xl = list(xn)
        poly = 0
        XProd = 1
        for i in range( len (a)):
            d = a[i]
            poly += d * XProd
            XProd = XProd * (x - xl[i])
        return (simplify(poly))


def DivDiffY(f,xi,yi,zi,M):
#function returns the lists of yn from interpolating f for fixed (x,z)
    y = symbols ('y')
    while (len(yi)) < 3:
        R = randint(0,M)
        if R not in yi:
            yi = np.append(yi,R)
    fn = np.array([int(f(xi,yi[i],zi)) for i in range (3)])
    yn = yi.copy()
    while True:
        m = len (yn)
        a = fn.copy()
        for k in range (1,m):
            a[k:m] = (a[k:m] - a[k-1:m-1])/(yn[k:m] - yn[0:m-k])
        while True:
            R = randint(0,M)
            if R not in yn:
                yn = np.append(yn,R)
                break
            else:
                continue
        fn = np.append(fn , f(xi,yn[-1],zi))
        if a[-1] ==0:
            return(yn)

def DivDiffX(f,xi,yi,zi,R):
#function returns the lists of xn from interpolating f for fixed (y,z)
    x = symbols ('x')
    ran = R
    xi = np.array([randint(1,ran)])
    while (len(xi)) < 3:
        R = randint(0,ran)
        if R not in xi:
            xi = np.append(xi,R)
    fn = np.array([int(f(xi[i],yi,zi)) for i in range (3)])
    xn = xi.copy()
    while True:
        m = len (xn)
        a = fn.copy()
        for k in range (1,m):
            a[k:m] = (a[k:m] - a[k-1:m-1])/(xn[k:m] - xn[0:m-k])
        while True:
            R = randint(0,ran)
            if R not in xn:
                xn = np.append(xn,R)
                break
            else:
                continue
        fn = np.append(fn , f(xn[-1],yi,zi))
        if a[-1] ==0:
            return(xn)
            
def Newton_3V(f, R): 
#combines DivDiffX, DivDiffY,PolyNewtonX and PolyNewtonY to return f(x,y,z)
    x ,y ,z = symbols ('x y z')
    ran = R
    xi  = np.array([randint(0,ran)])
    yi  = np.array([randint(0,ran)])
    zi  = np.array([randint(0,ran)])
    xn =  DivDiffX(f,xi,yi,zi,R)
    yn =  DivDiffY(f,xi,yi,zi,R)   
    XYM = []
    for i in range(len(xn)):
        XYM.append([PolyNewtonZ(f,xn[i],yn[j],R) for j in range(len(yn))])
    for k in range(len(XYM)):
        XM = [PolyNewtonY(yn,XYM[k]) for k in range(len(XYM))]
    return (PolyNewtonX( xn, XM))

   
def Plot_Failure(f):
    #produces plot for each #S probability of failure 
    S = [10+5*i for i in range (19)]
    l = []  #list of probabilities for each S
    for j in range (len (S)):
        i = 0
        c = 0 #tc = time counter, c = error counter , i = iteration counter
        N = 10**3
        for i in range (N):
            RECON1 = Newton_3V(f, S[j])
            RECON2 = Newton_3V(f, S[j])
            i += 1
            if RECON1!= RECON2:
                c += 1
        l.append(c/N)
    plt.plot(S,l,'x')
    plt.xlabel('#S')
    plt.ylabel('Probability of incorrect reconstruction')
    plt.show()

  
def Plot_Time(f): 
    #produces plot for each #S of time taken
    S = [20+1000*i for i in range (30)]  
    tl = [] #time list
    for j in range (len(S)):
        c = i = tc = 0 #tc = time counter, c = error counter , i = iteration counter
        for i in range (10**3):
            start = time.time()
            p = Newton_3V(f, S[j])
            elapsed = time.time()-start
            tc += elapsed
            i += 1
        tl.append(tc/10**3)
    plt.plot(S,tl,'x')
    plt.xlabel('#S')
    plt.ylabel('Time taken for reconstruction')
    plt.show()
   

