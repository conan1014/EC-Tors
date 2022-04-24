from sage.all import *
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import pickle
import primefac



def isPrime(a):
    '''fast prime checker'''
    return all(a % i for i in range(2, int(a)))

def sign(a):
    if a>=0:
        return '+'
    else:
        return ''
    
def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(3,P) if isPrime(i)]

def good_red(a,b,p):
    '''determines if p is a prime of good reduction for the EC'''
    if (4*a**3 + 27*b**2)%p == 0:
        return False
    else:
        return True

def ls(a, p):
    '''legendre symbol (a/p)'''
    r = a%p 
    if r==0:
        return 0
    if r==1:
        return 1
    factors = prime_factors_simple(r)
    final=1
    for x in factors:
        if x==2:
            final = final*(-1)**((p**2-1)/8)
        else:
            final = final*ls(p, x)*(-1)**(((p-1)/2)*((x-1)/2))
    return int(final)

def prime_factors_simple(n):
    '''prime factors in simple form'''
    if is_prime(n):
        return [n]
    primelst = primes_inc2(n)
    i = 0
    dividing = n
    factorlst = []
    while i < len(primelst):
        while dividing % primelst[i] == 0:
            z = primelst[i]
            factorlst.append(z)
            dividing = dividing/z
        i = i+1
    return factorlst

def prime_factors_simple_red(n):
    '''prime factors in simple form'''
    if is_prime(n):
        return [n]
    primelst = primes_inc2(n)
    i = 0
    dividing = n
    factorlst = []
    while i < len(primelst):
        while dividing % primelst[i] == 0:
            z = primelst[i]
            factorlst.append(z)
            dividing = dividing/z
        i = i+1
    return list(set(factorlst))

def primes_inc2(P):
    '''generate all primes < P'''  
    return [i for i in range(2,P) if isPrime(i)]

def valid_points_slow(a,b,x,p):
    y_2 = (x**3 + a*x + b)%p
    l_s = ls(y_2, p)
    if l_s >= 0:
        return 2*l_s
    else:
        return 0
    
def valid_points_slow_2(a,b,x,p):
    y_2 = (x**3 + a*x + b)%p
    QRS = [x**2%p for x in range(p)]
    if y_2%p==0:
        return 1
    elif y_2 in QRS:
        return 2
    else:
        return 0

def valid_points_quick(a,b,x,p):
    y_2 = (x**3 + a*x + b)%p
    if y_2%p==0:
        return 1
    elif QRs[p-1][y_2-1]==1:
        return 2
    else:
        return 0
    
def gp_size(a,b,p):
    '''Returns the number of order 3 points points in E[F_p]'''
    
    R = PolynomialRing(GF(p), 'x')
    x = R.gen()
    f = 3*x**4 + 6*a*x**2 + 12*b*x - a**2
    roots = f.roots()
    if len(roots)==0:
        return 1
    count = 1
    for n in roots:
        count += valid_points_quick(a,b,int(n[0]),p)
    return count


def gp_size_sage(a,b,p):
    '''Returns the number of order 3 points points in E[F_p] using sage'''
    
    F = GF(p)
    E = EllipticCurve(F, [a, b])
    P=E(0)
    return len(P.division_points(3))

def approx(val, frac, error):
    if frac-error<val<frac+error:
        return True 
    else:
        return False 
    
def type_prop(z, error):
    if approx(z[2], 27/48, error):
        return 1
    elif approx(z[2], 33/48, error):
        return 2
    elif approx(z[3], 2/3, error):
        return 3
    elif approx(z[4], 1/4, error):
        return 4
    elif approx(z[2], 0, error):
        return 5
    else:
        print("no classification", z)
        return

def gen_qrs(N):
    
    '''Generate QRs for all primes<N'''
    
    P = primes(N)
    Z = np.empty((1), int)[1:]

    
    for i in range(1,N+1):
        if i in P:
            Z = np.append(Z,[i],axis=0)
        else:
            Z = np.append(Z,[0],axis=0)
        
    QRs = np.zeros((N,N), int)
    
    for i in range(1,N+1):
        if Z[i-1]!=0:
            qrs = np.asarray(quadratic_residues(i))
            for x in qrs:
                QRs[i-1][x-1] = 1
        
    return QRs

tqr1 = time.time()
filename = 'QRs_10000.pk'
'''
QRs = gen_qrs(50000)
with open(filename, 'wb') as fi:
    # dump your data into the file
    pickle.dump(QRs, fi)
'''

with open(filename, 'rb') as fi:
    QRs = pickle.load(fi)
tqr2 = time.time()

print("QR generation time = ", tqr2-tqr1)


class EC_3tors:
    '''Class to calculated data and plots for 3-torsion group sizes'''
        
    def __init__(self, P, n):
        self.n = n
        self.P = P
        self.primes = primes(self.P)
        
    def single_curve_final(self, a, b):
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        data = [[p, gp_size(a,b,p)] for p in good_primes]
        no_1, no_3, no_9, tot = 0,0,0,0
        for arr in data:
            if arr[1]==1:
                no_1+=1
            elif arr[1]==3:
                no_3+=1
            else:
                no_9+=1
            tot+=1
        if no_1 + no_3 + no_9 != tot:
            print('ERROR')
            return
        return [a,b,no_1/tot,no_3/tot,no_9/tot]
    
    def single_curve_ls(self, a, b):
        no_1, no_3, no_9 = 0,0,0 #counters
        Z = np.empty((0,4), int)
        good_primes = [p for p in self.primes if good_red(a,b,p)]
        for p in good_primes:
            no = gp_size(a,b,p)
            if no==1:
                no_1+=1
            elif no==3:
                no_3+=1
            elif no==9:
                no_9+=1
            tot = no_1 + no_3 + no_9
            Z = np.append(Z, [[p, no_1/tot,no_3/tot,no_9/tot]], axis=0)
        return Z
    
    def plot_graph_single(self, a, b):
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        
        Z = self.single_curve_ls(a,b)
        Ps = Z[:,0]
        n1 = Z[:,1]
        n3 = Z[:,2]
        n9 = Z[:,3]
        print(n1[-1], n3[-1], n9[-1])
        plt.plot(Ps, n1, 'b-', label='#1')
        plt.plot(Ps, n3, 'r-',label='#3')
        plt.plot(Ps, n9, 'g-', label='#9')
        plt.ylabel('E[2] size ratio')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
        
    def mult_curve_data(self):
        n = self.n
        Z = np.empty((0,5), int)
        all_coeffs = itertools.product(range(-n,n+1), range(-n,n+1))
        coeffs = [[a,b] for a,b in all_coeffs if (4*a**3 + 27*b**2) != 0]
        
        for a,b in coeffs:
            z = self.single_curve_final(a,b)
            Z = np.append(Z, [[a,b,z[2],z[3],z[4]]], axis=0)
            if z[2] + z[3] + z[4]>1.0001 or z[2] + z[3] + z[4]<0.9999:
                print('error')
                return
        return Z
    
    def sort_data(self, error):
        Z = self.mult_curve_data()
        types = [type_prop(z,error) for z in Z]
        length = range(len(Z))
        type1 = [Z[i] for i in length if types[i]==1]
        type2 = [Z[i] for i in length if types[i]==2]
        type3 = [Z[i] for i in length if types[i]==3]
        type4 = [Z[i] for i in length if types[i]==4]
        type5 = [Z[i] for i in length if types[i]==5]
        types = [type1,type2, type3, type4, type5]
        return types
    
    def plot_mult(self, error):
        types = self.sort_data(error)
        type1 = np.asarray(types[0])
        type2 = np.asarray(types[1])
        type3 = np.asarray(types[2])
        type4 = np.asarray(types[3])
        type5 = np.asarray(types[4])
        plt.plot(type1[:,0], type1[:,1], 'y.')
        plt.plot(type2[:,0], type2[:,1], 'b.')
        plt.plot(type3[:,0], type3[:,1], 'k.')
        plt.plot(type4[:,0], type4[:,1], 'g.')
        plt.plot(type5[:,0], type5[:,1], 'm.')
        plt.grid()
        plt.xlabel('a')
        plt.ylabel('b')
        plt.show()
        
    def plot_scatter(self):
        Z = self.mult_curve_data()
        #print(Z[0:10])
        No_1 = Z[:,2]
        No_3 = Z[:,3]
        plt.scatter(No_1, No_3, marker='.')
        plt.title("Distribution scatter graph for P = {0}, n ={1}".format(self.P, self.n))
        plt.xlabel('#E[3] = 1 proportion')
        plt.ylabel('#E[3] = 3 proportion')
        
    def data_export(self):
        Z = self.mult_curve_data()
        filename_exp = 'Distrib_data_P={0}_n={1}.pk'.format(self.P, self.n)
        with open(filename_exp, 'wb') as fi:
            # dump your data into the file
            pickle.dump(Z, fi)
        
        
        
t0 = time.time()

'''
#GRAPHS TO PLOT IN PROJECT
E = EC_3tors(10000,20) 
E.plot_graph_single(2, -3)
E.plot_graph_single(2, 0)
E.plot_graph_single(0,3)
E.plot_graph_single(0,2)
E.plot_graph_single(0,-3)
E.plot_graph_single(-6,-8)
'''


E_mult = EC_3tors(1000,20) 
E_mult.plot_scatter()


t1 = time.time()



print('EC computation time = ', t1-t0)





