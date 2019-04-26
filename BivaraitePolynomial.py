"""
Bivariate Polynomial 

This file includes black box reconstruction algorithm in Newton_2Var(f,R) for a bivariate 
polynomial f and maximum random integer R

and the code for time plots and accuracy testing for reconstruction
of a bivariate polynomai in Plot_Time(f) and Test_Accuracy(f)

"""

import numpy as np 
from matplotlib import pyplot as plt
import time
from random import *
from fractions import Fraction
from sympy import *

x  = symbols ('x')
y  = symbols ('y')

def MatMul(X , Y):
    M = []
    result = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*Y)] for X_row in X]
    for r in result:
        M.append(r)
    return (M)
        

def Newton_2Var(f , R):
    #reconstructs a bivariate function f using th black box reconstruction formula
    ran = R
    xn = np.array([])
    yn = np.array([[0]])
    while (len(xn)) < 2:
        R = randint(0,ran)
        if R not in xn:
            xn = np.append(xn,R)      
    while (len(yn))<3:
        R = randint(0,ran)
        if R not in yn:
            yn = np.append(yn ,[[R]], axis =0)
    C = np.vstack((xn , np.zeros((2,2)))) 
    C = np.append(yn , C , axis =1)
    DS = 1
    n=0
    #makes matrix C for 2 random x , y values for f(x)
    while DS !=0 :
        DS = 0 #diagonal sum
        xn = C[0][1:]
        yn = C[1:,0]
        if n>0:
             #adds a new x and y values
            while True:
                Rx = randint(0,100)
                if Rx not in xn:
                    ax = np.array([[Rx]])
                    break
                else:
                    continue
            while True: 
                Ry = randint(0,100)
                if Ry not in yn:
                    ay = np.array([[Ry], [0]])
                    break
                else:
                    continue
            for m in range(len(C)-1):
                ax = np.vstack((ax,[0]))
                ay = np.append(ay , 0)
            C = np.append(C , ax , axis=1)
            C = np.vstack((C , [ay] ))
        n+=1
        #fills grid with function of x & y, in upper ∆
        for i in range ( 1, len(C) ):       
            for j in range (1, len(C)):
                if i+j > len(C[0]):
                    C[i][j] = 0
                else:
                    C[i][j] = f(C[0][j],C[i][0])
        #performs divided difference algorithm on C
        for k in range (1,len(C) ):
            for i in range ( len(C)-1, 0 , -1):       
                for j in range (len(C[0])-1 , 0 , -1):
                    if j <= (k) and i > (k):
                        C[i][j] = Fraction(int(C[i][j] - C[i-1][j]) , int(C[i][0] - C[i-k][0]))
                    if i <= (k) and j > (k):
                        C[i][j] = Fraction (int(C[i][j] - C[i][j-1]) , int(C[0][j] - C[0][j-k]))
                    if i > k and j > k:
                        C[i][j] = Fraction(int(C[i][j] + C[i-1][j-1] - C[i][j-1] - C[i-1][j]) , int((C[i][0] - C[i-k][0])*(C[0][j] - C[0][j-k])))
                    if i <= k and j <= k:
                        C[i][j] = C[i][j]
                    if i+j>len(C[0]):
                        C[i][j]=0
        for l in range ( 1, len(C)): 
            DS += C[l][len(C)-l]
    # Performs Y^tAX  
    X = [[1]]
    p = 1
    for i in range (len(xn)):
        p *= (x - xn[i])
        X.append([p])
    Y=[1]
    q =1
    for j in range (len(yn)):
        q *= (y - yn[j])
        Y.append(q)
    B =  [C[i][1:] for i in range (1,len (C))] #trim first row and column of x and y
    BX = MatMul ( B, X)
    final = MatMul([Y] , BX)
    final2 = [ expand(final[k][0]) for k in range (len(final))]
    return(final2)

def Plot_Time(f):
    #produces plot of average time taken for each #S
    S = [10+150*i for i in range (20)]
    t = []  #list of time for each S
    for j in range (len (S)):
        i  = 0 #iteration counter
        ts =time.time()
        for i in range (10**3):
            p = Newton_2Var(f, S[j])
            i += 1
        t.append((time.time() - ts)/10**3)
    plt.plot(S , t)
    plt.ylim([0,max(t)*1.2])
    plt.ylabel('Average time taken for reconstruction')
    plt.xlabel('Size of #S')
    plt.show()

def Test_Accuracy(f):
    #tests accuracy of reconstructing by perfomring infinte reconstructions of polynomial
    i = 0 #iteration counter
    c = 0 #failure counter
    while True:
        i+=1
        RECON1 = Newton_2Var(f, 10**3)
        RECON2 = Newton_2Var(f, 10**3)
        if RECON1 != RECON2:
            print('ERROR', RECON1 , RECON2 , i , c) #visual evidence of failure with counters
            c +=1
        if i%10**6==0:
            print (i , c) #visual of how many iteraitons have occured
    