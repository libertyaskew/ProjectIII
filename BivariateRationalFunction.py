"""
Bivariate Rational Function

This file includes black box reconstruction algorithm in Bivar_Thiele(f,M) for a bivariate 
rational fuction f and maximum random integer M
and the code for time plots and accuracy testing for reconstruction
of a bivariate ration fucntion in Plot_Time(f) and Test_Accuracy(f)

"""

from fractions import Fraction
from sympy import * 
import numpy as np
from random import *
import time
from matplotlib import pyplot as plt

x , y = symbols('x y')
def f(x,y):
    return( Fraction(( 4*x+y+2),(y+x+1)))
    
def MatMul(X , Y): #matrix multiplication as used in Mat2()
    result = np.zeros((len(X),len(Y[0])))
    for i in range(len(X)):
   # iterate through columns of Y
       for j in range(len(Y[0])):
       # iterate through rows of Y
           for k in range(len(Y)):
               result[i][j] += X[i][k] * Y[k][j]
    return(result)        
    
def RecDiffX(f,xi,yj,M):
    #constructs reciprocal divided difference for x with fixed y returns a list of nodes
    ran = M #M is maximum random integer
    yj = yj[0]
    while (len(xi)) < 3:
        while True:
            R = randint(1,ran)
            try: #ensures no poles
                f(R,yj)
                if R not in xi: #ensures no repeats
                    xi = np.append(xi,R)
                    break
            except ZeroDivisionError:
                continue
    xi = [int(xi[i]) for i in range (len(xi))]
    fn = np.array([f(xi[i],yj) for i in range (len(xi))])
    xn = xi.copy()
    c = []
    while True:
        m = len(xn)
        a = fn.copy()
        b = []
        for k in range (1,m):
            b.append(a[k:m-1].copy())
            if a[-1]==a[-2]:
                return( xn[:-1]) #returns list of nodes with the final one deleted as this is used in termination condition so will give a divde by zero error
            a[k:m] = [ Fraction((xn[i+k]-xn[i]),(a[i+k] - a[i+k-1])) for i in range (m-k)]
            if k>1:
                a[k:m]+=list(b[k-2])
        for k in range (2 , len(a)):
            a[k]-= a[k-2]
        while True:
            R = randint(1,ran)
            try: #ensures no poles
                f(R,yj)
                if R not in xn: #ensures no repeats
                    xn = np.append(xn,R)
                    fn = np.append(fn , f(xn[-1],yj))
                    break
            except ZeroDivisionError:
                continue
            
def RecDiffY(f,xi,yj, M):
    #constructs reciprocal divided difference for x with fixed y returns a list of nodes
    ran = M #M is maximum random integer
    xi = xi[0]
    while (len(yj)) < 4:
        while True:
            R = randint(1,ran)
            #print (x,R)
            try:
                f(xi,R)
                if R not in yj:
                    yj = np.append(yj,R)
                    break
            except ZeroDivisionError:
                continue
    yj = [int(yj[i]) for i in range (len(yj))]
    fn = np.array([f(xi,yj[i]) for i in range (len(yj))])
    yn = yj.copy()
    c = []
    while True:
        m = len(yn)
        a = fn.copy()
        b = []
        for k in range (1,m):
            b.append(a[k:m-1].copy())
            if a[-1]==a[-2]:
                return(yn[:-1])  #returns list of nodes with the final one deleted as this is used in termination condition so will give a divde by zero error
            a[k:m] = [ Fraction((yn[i+k]-yn[i]),(a[i+k] - a[i+k-1])) for i in range (m-k)]
            if k>1:
                a[k:m]+=list(b[k-2])
        for k in range (2 , len(a)):
            a[k]-= a[k-2]
        while True:
            R = randint(1,ran)
            try:
                f(xi,R)
                if R not in yn:
                    yn = np.append(yn,R)
                    fn = np.append(fn , f(xi,yn[-1]))
                    break
            except ZeroDivisionError:
                continue

def Bivar_Thiele(f,M):
#finds bivariate rational funciton reconstruction of f for maximum random integer M
    x  = symbols('x')
    y = symbols('y')
    ran = M #maximum random integer
    xi  = np.array([randint(0,ran)])
    yi  = [randint(0,ran)]
    xn = np.array(RecDiffX(f,xi,yi,M)) #calls list of nodes from univaraite reciprocal divded difference for x
    yn = RecDiffY(f,xi,yi,M)#calls list of nodes from univaraite reciprocal divded difference for y
    xn = np.append(0,xn)
    xn = np.transpose([(xn)])
    C = np.vstack((yn , np.zeros((len(xn)-1,len(yn)))))   
    C = np.append(xn , C , axis =1)
    C = np.array(C)
#finds reciprocal divded difference matrix
    k = 0   #iteration counter
    for i in range ( 1, len(C) ):       
        for j in range (1, len(C)):
            C = C.astype('object')
            C[i][j] = f(int(C[0][j]),int(C[i][0]))
    while k<5:
        k += 1 
        for i in range ( len(C)-1, 0 , -1):       
            for j in range (len(C[0])-1 , 0 , -1):
                if i+j<=k+1:
                    #print ('*',i,j,k)
                    C[i][j] = C[i][j]
                else:
                    if  j<= k :
                        C[i][j] = Fraction(int(C[i][0] - C[k-j+1][0]),(C[i][j]-C[k-j+1][j]))
                    if j>k:
                        C[i][j] = Fraction(int(C[0][j] - C[0][k]),(C[i][j]-C[i][k]))
    Y=[1]
    q =1
    xn = xn[1:]
    l = 0
    for l in range (len(yn)-1):    #construct Y vector 
        q *= (y - yn[l])
        Y.append(q)
    Y = np.array(Y)
    B =  np.array([C[i][1:] for i in range (1,len (C))])
    expr = 0
    for i in range(len(xn)-1,-1,-1): #works in reverse order as starts from bottom of continued fraction
        L = 0 #L is used to store A_j(y) for ea
        for j in range(len(Y)):
            L+= B[j][i]*Y[j] #multiplies Y vector by column of B to find A_j(y)
        if i == (len(xn)-1): 
            expr = L #last point in continued fraction
        else: #constructs 
            expr = (x - xn[i])/ expr
            expr+= L
        expr = simplify(expr) #sympy command to simplify expression
    return (expr) 

def Plot_Time(f):
    #produces plot of average time taken for each #S
    S = [10+2*i for i in range (20)]
    t = []  #list of time for each S
    for j in range (len (S)):
        i  = 0 #iteration counter
        ts =time.time()
        for i in range (10**3):
            p = Bivar_Thiele(f, S[j])
            i += 1
        t.append((time.time() - ts)/10**3)
    plt.plot(S , t)
    plt.ylim([0,max(t)*1.2])
    plt.ylabel('Average time taken for reconstruction')
    plt.xlabel('Size of #S')
    plt.show()

def Test_Accuracy(f):
    #tests accuracy of reconstructing by perfomring infinte reconstructions of bivaraite rational function f
    i = 0 #iteration counter
    c = 0 #failure counter
    while True:
        i+=1
        RECON1 = Bivar_Thiele(f, 50)
        RECON2 = Bivar_Thiele(f, 50)
        if RECON1 != RECON2:
            print('ERROR', RECON1 , RECON2 , i , c) #visual evidence of failure with counters
            c +=1
        if i%10**6==0:
            print (i , c) #visual of how many iteraitons have occured

Plot_Time(f)
Test_Accuracy(f)   