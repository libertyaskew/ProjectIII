"""
INTEGER RECONSTRUCTION

CRT_recon implements the integer reconstruction algorithm with 
the chinese remainder theorem as presented in the paper

Code for plots in paper is in Plot_Time and Plot_MaxInt and Plot_Accuracy
"""


import random
from numpy import *
from functools import reduce
from matplotlib import pyplot as plt
import time
import numpy as np
  
def ran_primes(n):
    #generates a random primes in range n
    noprimes = set(j for k in range(2, int(sqrt(n))+1) for j in range(k*2, n, k))
    primes =  [x for x in range(3, n) if x not in noprimes]
    return random.choice(primes)


def egcd(s, t):
    #returns r = xs + yt in form (r, x , y)
    if s == 0:
        return (t, 0, 1)
    else:
        g, y, x = egcd(t % s, s)
        return (g, x - (t // s) * y, y)

 
def mod_inv(a, m):
    #returns inverse a mod m
    g, x, y = egcd(a, m)
    if abs(g) != 1:
        raise Exception('mod inv does not exist')
    else:
        return x % m

def CRT(n , a):
    #returns chinese remainder of
    # x = a1 mod n1
    # x = a2 mod n2
    prod = np.prod(n)
    quo  = [ int(prod / n[i]) for i in range (len(n))]
    m_i  = [ ((mod_inv(quo[i] , n[i])) * quo[i] * a[i]) for i in range (len(n))]
    return( sum(m_i)%prod )
   
    
def CRT_recon(A , R):
#loop reconstructs A by taking its modulus from list of random primes n_i between
#1 & 100. Terminates when most recent 2 terms are equal. Works up to 13 int for A
#reconstructs A using CRT algorithm, returns estimation and number of iterations
    c = 2
    n_i = [ran_primes(R)]
    recon = [1]
    #generates a list of 2 unique primes
    if len(n_i) < 3:
        n_i.append(ran_primes(R))
        n_i = list(set(n_i))
    while True :
        c +=1
        rp = ran_primes(R)
        if rp not in n_i:
            n_i.append(rp)
        else:
            continue
        recon.append(int(CRT(n_i , [A%n_i[k] for k in range (len(n_i)) ] )))
        if recon[-1] == recon[-2]: #adding more recon equalities improves accuracy 
            return (recon[-1] , c)

def Plot_Time():
    #plots a scatter graph of average time taken for reconstruction
    #size of primes.
    tl = []
    cl = []
    mpl = []
    for i in range (100):
        A = random.randint(1, 10**6)
        MP = 200 + 13* i
        ts =time.time()
        for i in range(100):
            CR = CRT_recon(A,MP)
        tl.append(time.time() - ts)
        cl.append(CR[1])
        mpl.append(MP)
    plt.scatter(mpl , tl , marker = "x")
    plt.ylabel('Time taken for reconstruction')
    plt.xlabel('Size of maximum prime')
    plt.ylim([0,max(tl)*1.1])
    plt.savefig('TimeTaken.png')
    plt.show()


def Plot_Accuracy():
    #plot scatter graph showing difference between actual integer and constructed integer
    #and probability of failure for different primes
    ecl = []  #list of error count for each max prime
    mpl = []  #list of max primes
    rel = []  #list of relative error
    mpel = [] #list of max prime for each error
    for k in range (1000):
        MP = random.randint(200,1600)
        A = random.randint(1, 10**6)
        EC = RE = 0 #error counter & relative error
        for j in range (1000):
            CR = CRT_recon(A,MP) 
            if CR[0] != A:
                EC += 1
                RE += (CR[0] - A)/A
                rel.append(RE)
                mpel.append(MP)
            else:
                continue
        ecl.append(EC/1000)
        rel.append(RE / 1000)
        mpl.append(MP)
    plt.scatter((mpl) , (ecl) , marker = "x")
    plt.ylabel('Probability of incorrect reconstruction')
    plt.xlabel('Size of maximum prime')
    plt.ylim([0,max(ecl)*1.1])
    plt.savefig('ProbIncorrect.png')
    plt.show()
    plt.scatter((mpl) , (rel) , marker = "x")
    plt.ylim([0,max(rel)*1.2])
    plt.ylabel('Relative error')
    plt.xlabel('Size of maximum prime')
    plt.savefig('RelError.png')
    plt.show()


def Plot_MaxInt():
#plots maximum prime v max integer possible to reconstruct
    k = 0
    al = []
    pl = []
    for i in range (40):
        ma = 100
        ac = 0
        i +=1
        while True:
            try:
                CRT_recon(ma,100+10*i)
            except Exception:
                ac += ma
                al.append(ma)
                pl.append(100+10*i)
                break #failure when mod inverse 
            ma +=100  #increases size of integer by 100
    plt.plot(pl,al, 'x')
    plt.ylabel('Maximum Integer')
    plt.xlabel('Size of maximum prime')
    plt.savefig('MaxInt.png')
    

Plot_Time()
Plot_Accuracy()
Plot_MaxInt()



