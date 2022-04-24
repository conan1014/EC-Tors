import numpy as np
import matplotlib.pyplot as plt
from sage.all import *
from matplotlib import cm
import matplotlib.colors as clr
import itertools

'''small functions'''

def sign(a):
    if a>=0:
        return '+'
    else:
        return ''
 
def good_red(a,b,p):
    '''determines if p is a prime of good reduction for the EC'''
    if (4*a**3 + 27*b**2)%p == 0:
        return False
    else:
        return True
    
def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(3,P) if isPrime(i)]

def f(a,b,x):
    '''calculate f(x)'''
    return x**3 + a*x +b

def E_num2_primes(a,b,P,k):
    '''return all primes <P for which the legendre symbol for (-4*a-3*k**2/p)!=1,
    i.e. -4*a-3*k**2 a NQR mod p'''
    ps = primes(P)
    ps = [p for p in ps if ls(-4*a-3*k**2, p)!=1]
    return ps
    
def E_num4_primes(a,b,P,k):
    '''return all primes <P for which the legendre symbol for (-4*a-3*k**2/p)==1,
    i.e. -4*a-3*k**2 a QR mod p'''
    ps = primes(P)
    ps = [p for p in ps if ls(-4*a-3*k**2, p)==1]
    return ps

def isPrime(a):
    '''fast prime checker'''
    return all(a % i for i in range(2, int(a)))
   
def gens(a,b,p):
    E = EllipticCurve(GF(p), [a,b])
    return len(E.gens())
    
        
class Gp_struct:
    '''a group to generate structures of EC groups mod p for p<P'''
    
    def __init__(self, P, n):
       self.n = n
       self.P = P
       self.primes = primes(self.P)
      
    def single_ls(self, a, b):
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        cyclic, not_cyclic, triv = 0,0,0 #counters
        Z = np.empty((0,4), float)
        tot = 0
        for p in good_primes:
            E = EllipticCurve(GF(p), [a,b])
            no = len(E.gens())
            if no==1:
                cyclic=cyclic+1
            elif no==2:
                not_cyclic=not_cyclic+1
            else:
                triv=triv+1
            tot=tot+1
            Z = np.append(Z, [[p, cyclic/tot,not_cyclic/tot,triv/tot]], axis=0)      
        return Z
    
    def single_final(self, a, b):
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        data = [[p, gens(a,b,p)] for p in good_primes]
        cyclic, not_cyclic, triv, tot = 0,0,0,0
        for arr in data:
            if arr[1]==1:
                cyclic+=1
            elif arr[1]==2:
                not_cyclic+=1
            else:
                triv+=1
            tot+=1
        #assert(tot == cyclic+not_cyclic+triv)
        return [a,b,cyclic/tot,not_cyclic/tot,triv/tot]
    
    def plot_single(self, a, b):
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        
        Z = self.single_ls(a,b)
        Ps = Z[:,0]
        cyclic = Z[:,1]
        not_cyclic = Z[:,2]
        triv = Z[:,3]
        print(cyclic[-1], not_cyclic[-1], triv[-1])
        plt.plot(Ps, cyclic, 'b-', label='cyclic')
        plt.plot(Ps, not_cyclic, 'r-',label='not cyclic')
        plt.plot(Ps, triv, 'g-', label='trivial')
        plt.ylabel('cyclicity proportions')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
        return
    
    def mult_data(self):
        n = self.n
        Z = np.empty((0,5), int)
        all_coeffs = itertools.product(range(-n,n+1), range(-n,n+1))
        coeffs = [[a,b] for a,b in all_coeffs if (4*a**3 + 27*b**2) != 0]
        
        for a,b in coeffs:
            z = self.single_final(a,b)
            Z = np.append(Z, [[a,b,z[2],z[3],z[4]]], axis=0)

        return Z
    
    def plot_scatter(self):
        Z = self.mult_data()
        cyc = Z[:,2]
        not_cyc = Z[:,3]
        plt.scatter(cyc, not_cyc, marker='.')
        plt.title("Cyclicity scatter graph for P = {0}, n ={1}".format(self.P, self.n))
        plt.xlabel('Cyclic proportion')
        plt.ylabel('Not cyclic proportion')
        return self
  

cyc = Gp_struct(5000, 5)

cyc.plot_scatter()


#print(cyc.plot_single(0,7))


      































