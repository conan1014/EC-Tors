import numpy as np
from sage.all import *
import matplotlib.pyplot as plt



def fact_modp(coeff, p):
    '''factorises polynomial into irreducibles mod p'''
    
    R = PolynomialRing(GF(p), 'x')
    x = R.gen()
    g = 0
    for i in range(len(coeff)):
        g += coeff[i]*x**(i)
    return list(g.factor())

def roots_modp(coeff, p):
    factors = fact_modp(coeff, p)
    roots = sum([factors[i][1] for i in range(len(factors))])
    return roots

def fact_Q(coeff):
    '''factorises polynomial into irreducibles mod p'''
    
    R = PolynomialRing(RationalField(), 'x')
    x = R.gen()
    g = 0
    for i in range(len(coeff)):
        g += coeff[i]*x**(i)
    return list(g.factor())

def roots_Q(coeff):
    factors = fact_Q(coeff)
    lin_factors = [factors[i] for i in range(len(factors)) if factors[i][0].degree()==1]
    roots = sum([lin_factors[i][1] for i in range(len(lin_factors))])
    return roots

def quart_roots(a,b):
    coeff = [3, 0, 6*a, 12*b, -a**2]
    #coeff = [-a**2, 12*b, 6*a, 0, 3]
    return roots_Q(coeff)

def roots_data(n):
    Z = np.empty((0,3), int)
    all_coeffs = itertools.product(range(-n,n+1), range(-n,n+1))
    coeffs = [[a,b] for a,b in all_coeffs if (4*a**3 + 27*b**2) != 0]
    
    for a,b in coeffs:
        z = quart_roots(a,b)
        if z not in [0,1,2]:
            print(z)
        Z = np.append(Z, [[a,b,z]], axis=0)
    return Z

def plot_struct(n):
    Z = roots_data(n)
    n0 = np.asarray([z for z in Z if z[2]==0])
    n1 = np.asarray([z for z in Z if z[2]==1])
    n2 = np.asarray([z for z in Z if z[2]==2])
    n3 = np.asarray([z for z in Z if z[2]==4])
    print(n3)
    plt.plot(n0[:,0], n0[:,1], 'y.')
    plt.plot(n1[:,0], n1[:,1], 'b.')
    plt.plot(n2[:,0], n2[:,1], 'm.')
    plt.grid()
    plt.xlabel('a')
    plt.ylabel('b')
    plt.show()
    

plot_struct(30)

































