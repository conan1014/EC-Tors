from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns

'''
THIS CLASS USES SAGE
'''

def sign(a):
    if a>=0:
        return '+'
    else:
        return ''
    
def isPrime(a):
    '''fast prime checker'''
    return all(a % i for i in range(2, int(np.sqrt(a))+1))

def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(3,P) if isPrime(i)]

def good_red(a,b,p):
    '''determines if p is a prime of good reduction for the EC'''
    if (4*a**3 + 27*b**2)%p == 0:
        return False
    else:
        return True
    
    
class ST:
    
    def __init__(self, P):
        self.P = P
        self.primes = primes(self.P)
        
    def mp(self, a, b):
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        Z = np.empty((0,2), int)
        for p in good_primes:
            E = EllipticCurve(GF(p), [a,b])
            Z = np.append(Z, [[p, E.order()]], axis=0)   
        return Z
    
    def theta(self, a, b):
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        Z = np.empty((0,3), int)
        for p in good_primes:
            E = EllipticCurve(GF(p), [a,b])
            mp = E.order()
            theta_ = math.acos((mp-p-1)/(2*np.sqrt(p)))
            Z = np.append(Z, [[p, mp, theta_]], axis=0)  
        return Z
    
    def plot_func(self, xmin, xmax):
        X = np.linspace(xmin, xmax, 1000)
        Y = [2/pi*math.sin(x)**2 for x in X]
        plt.plot(X,Y,label='A')
        
    def plot_hasse1(self, a, b):
        mps = self.mp(a,b)
        X = np.linspace(0, self.P, 1000)
        core = [y+1 for y in X]
        ymin = [y + 2*y**0.5 for y in core]
        ymax = [y - 2*y**0.5 for y in core]
        plt.plot(mps[:,0], mps[:,1], '.', label='#E(F_p)')
        plt.plot(core, ymin, 'r-')
        plt.plot(core, ymax, 'r-', label='Error bound')
        plt.xlabel('p')
        plt.ylabel('#E(F_p)')
        plt.legend()
        
    def plot_hasse2(self, a, b):
        mps = self.mp(a,b)
        X = np.linspace(0, self.P, 1000)
        core = [y+1 for y in X]
        ymin = [y + 2*y**0.5 for y in core]
        ymax = [y - 2*y**0.5 for y in core]
        mps_norm = [x[1]-x[0]-1 for x in mps]
        sns.distplot(mps_norm)
        #plt.plot(core, ymin, 'r-')
        #plt.plot(core, ymax, 'r-', label='Error bound')
        #plt.xlabel('p')
        #plt.ylabel('#E(F_p)')
        plt.legend()
        
        
    def plot_theta(self, a, b):
        Z = self.theta(a,b)
        thetas = Z[:,2]
        xmin, xmax = min(thetas),max(thetas)
        self.plot_func(xmin, xmax)
        sns.distplot(thetas, bins = 20, label='Theta_p distribution')
        plt.xlabel('x')
        plt.ylabel('y/density')
        plt.legend()
      
    def integrand(t):
        return (math.sin(t))**2
    
    def integral(self, alpha,beta):
        '''returns integral calculation'''
        
        result = quad(s_t_comparison.integrand,alpha,beta)
        return (2/math.pi)*result[0]
    
    def frac(self, alpha, beta, a, b):
        
        Z = self.theta(a, b)
        Z_red = [x for x in Z if alpha<= x[2] <=beta]
        return len(Z_red)/len(Z)

st = ST(1000)
st.plot_theta(1,1)
plt.show()
st.plot_hasse1(1,1)     
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            