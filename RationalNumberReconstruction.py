
"""
RATIONAL NUMBER RECONSTRUCTION

CRT_rat_recon implements the rational integer reconstruction algorithm combined with
the chinese remainder theorem as presented in the paper

Algorithm for producing efficency and maximum possible reconstruction plots is also included in
Plot_Time and Plot_MaxInt
"""

import numpy as np
from random import *
from functools import reduce
import time
from matplotlib import pyplot as plt

def egcd(s, t):
    #returns r = xs + yt in form (r, x , y)
    if s == 0:
        return (t, 0, 1)
    else:
        g, y, x = egcd(t % s, s)
        return (g, x - (t // s) * y, y)

def ran_primes(m,n):
    #generates a random primes in range n
    noprimes = set(j for k in range(2, int(np.sqrt(n))+1) for j in range(k*2, n, k))
    primes =  [x for x in range(m, n) if x not in noprimes]
    return choice(primes)
 
def mod_inv(a, m):
    #returns inverse a mod m
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def ZP2Q(  c,m):
    #returns rational number from integer in finite field mod m
    u = [1,0,m]
    v = [0,1,c]
    
    while np.sqrt(m/2)<=v[2]:
        q = int(u[2]/v[2])
        r = [(u[i] - q*v[i]) for i in range (3)]
        u = v
        v = r
    if abs(v[1])>=np.sqrt(m/2):
            raise Exception ('Error')
    return (v[2],(v[1]))

def Q2ZP( q , r , p):
    #returns rational number (q/r) as integer mod p
    return ((q * mod_inv(r,p))%p)


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
        #print(n_i , 'n-i')
        if recon[-1] == recon[-2]== recon[-3]:
            return (recon[-1] , c)

        
def CRT_rat_recon(MR , MP):
#reconstructs a rational number from random s and t by taking multiple evaluations in finite fields
#then implements the reconstruction algorithm presented in paper.
    c = 0 #iteration counter
    q  = ran_primes(3,MR) 
    r  = ran_primes(3,MR)
    if q == r:
        q = r = 1  
    pB = q
    while q%pB == 0 and r%pB == 0:
        pB = ran_primes(50,MP)
    p = [pB] #masterlist of primes to ensure no repeats
    ZB = Q2ZP(q,r,pB)
    while True:
        try:
            QB = ZP2Q(ZB,pB)
            break
        except Exception:
            QB = 0
            break
    while True: 
        c +=1
        while True:
            pA = ran_primes(50,MP)
            if pA not in p and q%pA != 0 and r%pA != 0:
                p.append(pA)
                break
            else:
                continue
        ZA = Q2ZP(q,r,pA)
        ZB = CRT([pA,pB],[ZA,ZB])
        pB = np.prod([pA,pB])
        while True:
            try:
                QA = ZP2Q(ZB, pB)
                break
            except Exception:
                QA = 0 #exception raised if prime is not big enough. by setting QA=0 effectively skips this iteration and waits until prime isbigger on next iteration. 
                break
        if QA == QB and QA!=0:
            if QA[0]!= q and QA[1]!= r:
                print(QA,QB, q ,r)
                break
            else:
                return(c)
                break
        else:
            QB = QA
            continue 

        
def Plot_MaxInt():
#produces plot of maximum size of s and t which can be reconstructed for a maximum prime
    i = 0
    k = 0
    pl = []
    rl = []
    for i in range (200):
        mr = 30
        while True:
            try:
                CRT_rat_recon(mr,100+10*i)
            except Exception: #continues until Exception error is raised from mod_inv function as integers are so large the computer registers them as 0.
                pl.append(100+10*i)
                rl.append(mr)
                break 
            mr+=10  #increases s and t by 10
        i +=1 
    print (pl,rl)
    plt.plot(pl,rl)
    plt.ylabel('Maximum values of q , r')
    plt.xlabel('Size of maximum prime')
    plt.savefig('RatRecon_MaxQR.png')
    plt.show()
    
    
def Plot_Time():
#produces plot for average time taken to recconstruct rational numbers
    i = 0
    R = 500 #fixed integer
    pl = []
    tl= []
    mp = 100
    while True:
        i = 0
        tc = 0
        try:
            for i in range(100):
                ts =time.time()
                CRT_rat_recon(R,mp)
                tc += (time.time() - ts)
                i+=1
            tl.append(tc)
            pl.append(mp)
        except Exception: 
            break #continues increasing prime and timing until failure
        mp+=5
        i +=1 
    plt.plot(pl,tl)
    plt.ylabel('Time taken for reconstructions')
    plt.xlabel('Size of maximum prime')
    plt.savefig('RatRecon_Time.png')
    plt.show()
    
ScatRatRecon1()
ScatRatRecon2()    


    




