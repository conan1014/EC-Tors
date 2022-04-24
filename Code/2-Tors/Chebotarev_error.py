import numpy as np
import matplotlib.pyplot as plt
from sage.all import *
from matplotlib import cm
import matplotlib.colors as clr


def sign(a):
    if a>=0:
        return '+'
    else:
        return ''


def isPrime(a):
    '''fast prime checker'''
    return all(a % i for i in range(2, int(a)))


def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(3,P) if isPrime(i)]


def roots_modp(a,b,p):
    '''determines number of second order points to the EC mod p using a fast algorithm'''
    R = PolynomialRing(GF(p), 'x')
    x = R.gen()
    f = x**p -x
    g = x**3 + a*x + b
    div = f.gcd(g)
    return div.degree()


def good_red(a,b,p):
    '''determines if p is a prime of good reduction for the EC'''
    if (4*a**3 + 27*b**2)%p == 0:
        return False
    else:
        return True


def no_elts_prop(a, b, P):
    '''run through all good primes to find size of E[2] for each reduction, 
    culmulative count given'''
    all_primes = primes(P)
    good_primes = [p for p in all_primes if good_red(a,b,p)]
    no_1 = 0 #counter
    no_2 = 0 #counter
    no_4 = 0 #counter
    Z = np.empty((0,4), int) #primes covered, #(size 2), #(size 4) 
    for p in good_primes:
        no = roots_modp(a,b,p)+1
        if no==1:
            no_1+=1
        elif no==2:
            no_2+=1
        elif no==4:
            no_4+=1
        tot = no_1 + no_2 + no_4
        Z = np.append(Z, [[p, no_1/tot, no_2/tot, no_4/tot]], axis=0) 
    return Z


def error(x, D, lnGM):
    lnx = np.log(x)
    final = (D**0.5*lnx*(lnx + lnGM))/x**0.5
    #final = lnx**2/x**0.5
    return final


def plot_error(n, D, G, M):
    X = np.linspace(20, n, n)
    lnGM = np.log(G*M)
    Y = [error(x, D, lnGM)-0.5 for x in X]
    print(Y[-1])
    plt.plot(X,Y)
    


def plot_graph_prop(a, b, P, D, G, M):
    if (4*a**3 + 27*b**2) == 0:
        print("Singular")
        return
    
    Z = no_elts_prop(a,b,P)
    Ps = Z[:,0]
    n1 = Z[:,1]
    n2 = Z[:,2]
    n4 = Z[:,3]
    print(n1[-1], n2[-1], n4[-1])
    #plt.plot(Ps, n1, 'b-', label='E[2]=1')
    #plt.plot(Ps, n2, 'r-',label='E[2]=2')
    #plt.plot(Ps, n4, 'g-', label='E[2]=4')
    plt.ylabel('Chebotarev density error term')
    plt.xlabel('p')
    plt.legend()
    #plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
    plot_error(P, D, G, M)
    plt.show()


plot_graph_prop(1, 5, 10000, 3, 6, 31)

































