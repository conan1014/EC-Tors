import numpy as np
import matplotlib.pyplot as plt
from sage.all import *
from matplotlib import cm
import matplotlib.colors as clr

'''small functions'''

def sign(a):
    if a>=0:
        return '+'
    else:
        return ''
    
def factors(n):
    '''return all factors of n'''
    ls = [n, -n]
    for x in range(1,int(n**0.5)+3):
        if n%x==0:
            ls.append(x)
            ls.append(-x)
    return ls

def no_int_roots(a,b):
    '''returns number of integer roots to f(x) = x^3+ax+b'''
    if b==0:
        factors1 = factors(abs(a))
        count=1
        factors1 = list(set(factors1))
        for x in factors1:
            if x**2 + a == 0:
                count+=1
        return count
    factors1 = factors(abs(b))
    factors1 = list(set(factors1))
    count=0
    for x in factors1:
        if x**3 + a*x + b == 0:
            count+=1
    return count 
    
def good_red(a,b,p):
    '''determines if p is a prime of good reduction for the EC'''
    if (4*a**3 + 27*b**2)%p == 0:
        return False
    else:
        return True
    
def is_prime(m):   
    if m==1:
        return False
    f = 2
    while True and f*f <= m:
        if m%f == 0:
            return False
        f = f+1
    return True
    
def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(3,P) if isPrime(i)]

def primes_inc2(P):
    '''generate all primes < P'''  
    return [i for i in range(2,P) if isPrime(i)]

def gen_root_data(n):
    '''generate number of rational roots to cubics for all combinations of coefficients 
    less than n'''
    Z = np.empty((0,3), int)
    for a in range(-n,n):
        for b in range(-n,n):
            if (4*a**3 + 27*b**2) != 0:
                Z = np.append(Z, [[a, b, no_int_roots(a,b)]], axis=0) 
    return Z

def gen_root_graph(n):
    
    '''plots rational roots to cubics with coefficients less than n'''
    Z = gen_root_data(n)
    types = [[],[],[],[]]
    for x in Z:
        types[x[2]].append(list(x))
    no_0 = np.asarray(types[0])
    no_1 = np.asarray(types[1])
    no_2 = np.asarray(types[2])
    no_3 = np.asarray(types[3])
    plt.plot(no_0[:,0], no_0[:,1], 'y.')
    plt.plot(no_1[:,0], no_1[:,1], 'b.')
  #  plt.plot(no_2[:,0], no_2[:,1], 'b.')
    plt.plot(no_3[:,0], no_3[:,1], 'r.')
    plt.xlabel('a')
    plt.ylabel('b')
    
    
    x = np.linspace(-n,n,100)
    Y_1 = [x*i+i**3 for i in range(10)]
    Y_2 = [-y for y in Y_1]
    Y = Y_1 + Y_2

    for y in Y:
        x_clip=[]
        y_clip=[]
        for i in range(len(y)):
            if abs(y[i]) <= n:
                x_clip.append(x[i])
                y_clip.append(y[i])
        plt.plot(x_clip,y_clip,'-k')
    
    
    plt.grid()
    plt.plot()
    
def f(a,b,x):
    '''calculate f(x)'''
    return x**3 + a*x +b

def prime_factors(n):
    '''prime factors with multiplicity form'''
    if is_prime(n):
        return [[n,1]]
    primelst = primes_inc2(n)
    i = 0
    dividing = n
    factorlst = []
    while i < len(primelst):
        while dividing%primelst[i] == 0:
            z = primelst[i]
            factorlst.append(z)
            dividing = dividing/z
        i = i+1
        
    flat_factors = []
    [flat_factors.append(x) for x in factorlst if x not in flat_factors]
    factorlst = np.array(factorlst)
    factors = [[x, np.count_nonzero(factorlst==x)] for x in flat_factors]
    return factors  

def prime_factors_simple(n):
    '''prime factors in simple form'''
    if is_prime(n):
        return [n]
    primelst = primes_inc2(n)
    i = 0
    dividing = n
    factorlst = []
    while i < len(primelst):
        while dividing%primelst[i] == 0:
            z = primelst[i]
            factorlst.append(z)
            dividing = dividing/z
        i = i+1
    return factorlst

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
   
def roots_modp(a,b,p):
    '''determines number of second order points to the EC mod p using a fast algorithm'''
    R = PolynomialRing(GF(p), 'x')
    x = R.gen()
    f = x**p -x
    g = x**3 + a*x + b
    div = f.gcd(g)
    return div.degree()

def good_ratio(x,y):
    if x==0 and y==0:
        return 0.5
    elif y==0:
        return 1
    elif x/y>1:
        return x/y
    else:
        return x/y
    

    
class EC_modp:
    """EC_modp: class of ECs mod p prime (y^2 = x^3 + ax + b)"""
    
    def __init__(self, a, b, p):
      self._c = [a,b] #coeffs
      self.p = p #mod p
      
    def ls(self, a):
        '''legendre symbol'''
        ls = pow(int(a), (self.p-1)//2, self.p)
        if ls == self.p:
            return -1
        else:
            return ls 
      
    def genpoints(self):
        ''' Generate all points in EC mod p '''
        a, b= self._c[0], self._c[1]
        x_poss = np.arange(self.p)
 
        QRs = [x**2%self.p for x in x_poss]
        CFp = np.array([[0,1,0]])
        
        for x in x_poss:
            
            y = (x**3 + a*x + b)%self.p
            legendre_symbol = self.ls(y)

            if legendre_symbol == 1:
                CFp = np.append(CFp, [[x, QRs.index(y), 1]], axis=0)
                CFp = np.append(CFp, [[x, (-QRs.index(y))%self.p, 1]], axis=0)
            
            if y == 0:
                CFp = np.append(CFp, [[x, y, 1]], axis=0)
         
        return(CFp)
    
    def gen_table(self):
        '''Generate table with all possible points mod p, error given when points not possible'''
        
        a, b = self._c[0], self._c[1]
        x_poss = np.arange(self.p) #all possible values for x
 
        QRs = [x**2%self.p for x in x_poss] #list of all QRs
        CFp_table = np.array([[0,1,1,0]]) 
        
        for x in x_poss:
            y = (x**3 + a*x + b)%self.p #possible y^2 value for given x
            legendre_symbol = self.ls(y) #legendre symbol of y^2
            
            if legendre_symbol == 1:
                '''append (x, plus or minus y) to CFp_table if y^2 is a QR'''
                CFp_table = np.append(CFp_table, [[x, QRs.index(y), y,1]], axis=0) 
                CFp_table = np.append(CFp_table, [[x, -QRs.index(y), y,1]], axis=0)
            
            elif y == 0:
                '''append to (x,0) if y^2=0'''
                CFp_table = np.append(CFp_table, [[x, y, y, 1]], axis=0)
                
            else:
                '''append error point for y^2 being a non quadratic residue (no root)'''
                CFp_table = np.append(CFp_table, [[x, EOFError, y, 1]], axis=0)
        
        return(CFp_table)
    

    
    def add(self, p1, p2):
        '''calculate p1 + p2'''
        a, b = self._c[0], self._c[1]
        x1, y1, z1 = p1[0], p1[1], p1[2]
        x2, y2, z2 = p2[0], p2[1], p2[2]
        if z1 == 0:
            x3 = x2
            y3 = y2
            z3 = z2
        elif z2 ==0:
            x3 = x1
            y3 = y1
            z3 = z1
        elif y1 == 0 and y2==0 and x1==x2:
            x3 = 0
            y3 = 1
            z3 = 0
        elif y1==y2 and x1==x2 and y1 != 0:
            m = (3*x1**2 + a)*pow(int(2*y1), -1, self.p)
            v = -m*x1 + y1
            x3 = int((m**2 -2*x1)%self.p)
            y3 = int((-(m*x3 + v))%self.p)
            z3 = 1
        elif x1==x2 and y1 == int((-y2)%self.p):
            x3 = 0
            y3 = 1
            z3 = 0
        else:
            m = (y2-y1)*pow(int(x2-x1), -1, self.p)
            v = y1 - m*x1
            x3 = int((m**2 -x1-x2)%self.p)
            y3 = int((-(m*x3 + v))%self.p)
            z3 = 1
        return [x3,y3, z3]
    
    def negative(self, p1):
        x1, y1, z1 = p1[0], p1[1], p1[2]
        return [x1, int((-y1)%self.p), 1]
        
    def order(self, x):
        count = 1
        temp = x
        while temp[2] != 0:
            temp = self.add(temp, x)
            count+=1
        return count
    
    def elts_ord_2(self):
        a, b = self._c[0], self._c[1]
        x_poss = np.arange(self.p) #all possible values for x
        x_elts = [x for x in x_poss if f(a,b,x)%self.p==0] 
        elts = [[x,0,1] for x in x_elts]
        elts.append([0,1,0])
        return elts
        


 
    
class E_2:
    ''''a group to generate number of elements of order 2 mod p, for all p 
    prime p<P'''
    
    def __init__(self, P):
        self.P = P
    
    def no_elts_ls(self, a, b):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given'''
        all_primes = primes(self.P)
        good_primes = [p for p in all_primes if good_red(a,b,p)]
        no_1 = 0 #counter
        no_2 = 0 #counter
        no_4 = 0 #counter
        Z = np.empty((0,4), int) #primes covered, #(size 2), #(size 4) 
        for p in good_primes:
            ec = EC_modp(a,b,p)
            ord_2 = ec.elts_ord_2()
            no = len(ord_2) #size of group E[2]
            if no==1:
                no_1+=1
            elif no==2:
                no_2+=1
            elif no==4:
                no_4+=1
            Z = np.append(Z, [[p, no_1, no_2, no_4]], axis=0) 
        return Z
    
    def no_elts_ls_fast(self, a, b):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given'''
        all_primes = primes(self.P)
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
            Z = np.append(Z, [[p, no_1, no_2, no_4]], axis=0) 
        return Z
    
    def no_elts_prop(self, a, b):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given'''
        all_primes = primes(self.P)
        good_primes = [p for p in all_primes if good_red(a,b,p)]
        print(len(good_primes))
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
    
    def no_elts_ratios(self, a, b):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given - DO FOR ONLY GRAPH TYPE 1'''
        if no_int_roots(a,b)!=0:
            print('incorrect graph type')
            return 
        all_primes = primes(self.P)
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
            n1_n2 = good_ratio(no_1, no_2)
            n4_n1 = good_ratio(no_4, no_1)
            n4_n2 = good_ratio(no_4, no_2)
            Z = np.append(Z, [[p, n1_n2, n4_n1, n4_n2]], axis=0) 
        print(no_1, no_2, no_4)
        return Z
    
    def no_elts_ratios_normalised(self, a, b):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given - DO FOR ONLY GRAPH TYPE 1'''
        if no_int_roots(a,b)!=0:
            print('incorrect graph type')
            return 
        all_primes = primes(self.P)
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
            n1_n2 = good_ratio(no_1, no_2)
            n4_n1 = good_ratio(no_4, no_1)
            n4_n2 = good_ratio(no_4, no_2)
            Z = np.append(Z, [[p, n1_n2-2/3, n4_n1-1/2, n4_n2-1/3]], axis=0) 
        return Z
    
    def no_elts_ls_choose_primes(self, a, b, all_primes):
        '''run through all good primes to find size of E[2] for each reduction, 
        culmulative count given, choose primes'''
        good_primes = [p for p in all_primes if good_red(a,b,p)]
        no_2 = 0 #counter
        no_4 = 0 #counter
        no_1 = 0 #counter
        Z = np.empty((0,4), int) #primes covered, #(size 2), #(size 4) 
        for p in good_primes:
            ec = EC_modp(a,b,p)
            ord_2 = ec.elts_ord_2()
            no = len(ord_2) #size of group E[2]
            if no==2:
                no_2+=1
            elif no==4:
                no_4+=1
            elif no==1:
                no_1+=1
            Z = np.append(Z, [[p, no_1, no_2, no_4]], axis=0) 
        return Z
    
    def plot_graph_no_elts_choose_primes(self, a, b, all_primes):
        '''plot #E[2] for increasing p, choose primes'''
        
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        
        Z = E_2.no_elts_ls_choose_primes(self, a, b, all_primes) 
        Ps = Z[:,0]
        no_1 = Z[:,1]
        no_2 = Z[:,2]
        no_4 = Z[:,3]

        plt.step(Ps, no_1, 'b-', label='#E[2]=1')
        plt.step(Ps, no_2, 'r-',label='#E[2]=2')
        plt.step(Ps, no_4, 'g-', label='#E[2]=4')
        plt.ylabel('#E[2] count')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
    
    def plot_graph_no_elts(self, a, b):
        '''plot #E[2] for increasing p'''
        
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        
        Z = E_2.no_elts_ls(self, a, b) 
        Ps = Z[:,0]
        no_1 = Z[:,1]
        no_2 = Z[:,2]
        no_4 = Z[:,3]
        tot = no_1[-1] + no_2[-1] + no_4[-1]
        print(no_1[-1]/tot,no_2[-1]/tot,no_4[-1]/tot)

        plt.step(Ps, no_1, 'b-', label='#E[2]=1')
        plt.step(Ps, no_2, 'r-',label='#E[2]=2')
        plt.step(Ps, no_4, 'g-', label='#E[2]=4')
        plt.ylabel('Culmulative #E[2] count')
        
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
    
    def plot_graph_ratios(self, a, b):
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        if no_int_roots(a,b)!=0:
            print('incorrect graph type')
            return 
        
        Z = E_2.no_elts_ratios(self, a,b)
        Ps = Z[:,0]
        n1_n2 = Z[:,1]
        n4_n1 = Z[:,2]
        n4_n2 = Z[:,3]
        print(n1_n2[-1], n4_n1[-1], n4_n2[-1])
        plt.plot(Ps, n1_n2, 'b-', label='#1/#2')
        plt.plot(Ps, n4_n1, 'r-',label='#4/#1')
        plt.plot(Ps, n4_n2, 'g-', label='#4/#2')
        plt.ylabel('E[2] size ratio')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
        
    def plot_graph_props(self, a, b):
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return

        
        Z = E_2.no_elts_prop(self, a,b)
        Ps = Z[:,0]
        n1 = Z[:,1]
        n2 = Z[:,2]
        n4 = Z[:,3]
        print(n1[-1], n2[-1], n4[-1])
        plt.plot(Ps, n1, 'b-', label='E[2]=1')
        plt.plot(Ps, n2, 'r-',label='E[2]=2')
        plt.plot(Ps, n4, 'g-', label='E[2]=4')
        plt.ylabel('E[2] size proportions')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
        
    def plot_graph_ratios_normalised(self, a, b):
        if (4*a**3 + 27*b**2) == 0:
            print("Singular")
            return
        if no_int_roots(a,b)!=0:
            print('incorrect graph type')
            return 
        
        Z = E_2.no_elts_ratios_normalised(self, a,b)
        Ps = Z[:,0]
        n1_n2 = Z[:,1]
        n4_n1 = Z[:,2]
        n4_n2 = Z[:,3]
        print(n1_n2[-1], n4_n1[-1], n4_n2[-1])
        plt.plot(Ps, n1_n2, 'b-', label='#1/#2')
        plt.plot(Ps, n4_n1, 'r-',label='#4/#1')
        plt.plot(Ps, n4_n2, 'g-', label='#4/#2')
        plt.ylabel('E[2] size ratio')
        plt.xlabel('p')
        plt.legend()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(a), a, sign(b), b))
        plt.show()
        
    def ratio_check(self, n, error):
        Z = np.empty((0,4), int)
        for a in range(-n,n):
            for b in range(-n,n):
                if no_int_roots(a,b)==0:
                    if (4*a**3 + 27*b**2) != 0:
                        Z = E_2.no_elts_ratios_normalised(self, a, b)
                        print(Z[:,1][-1], Z[:,2][-1], Z[:,2][-1])
                        ls = np.array([Z[:,1][-1], Z[:,2][-1], Z[:,3][-1]])
                        ls = np.absolute(ls)
                        if np.amax(ls)>error:
                            print("y^2 = x^3 {0} {1}x {2} {3} gives significant error".format(sign(a), a, sign(b), b))
                            E_2.plot_graph_ratios_normalised(self, a, b)
                    
    
    def graph_type(self, a,b):
        
        '''determine graph type of input EC for increasing p'''
        
        Z = E_2.no_elts_ls(self, a, b)
        final = Z[-1]
        no_1, no_2, no_4 = final[1], final[2], final[3]
        if no_1>0 and no_2>0 and no_4>0:
            return 1
        elif no_1>0 and no_4>0 and no_2==0:
            return 1
        elif no_4>0 and no_1==0 and no_2==0:
            return 3
        elif no_4>0 and no_2>0 and no_1==0:
            return 4
        else:
            print("unrecognised graph type, try increasing P")
            print("unrecognised graph: a=", a,"b= ", b)
            return 6
        
    def graph_type_fast(self, a,b):
        '''determine graph type of input EC for increasing p'''
        
        Z = E_2.no_elts_ls_fast(self, a, b)
        final = Z[-1]
        no_1, no_2, no_4 = final[1], final[2], final[3]
        if no_1>0 and no_2>0 and no_4>0:
            return 1
        elif no_1>0 and no_4>0 and no_2==0:
            return 1
        elif no_4>0 and no_1==0 and no_2==0:
            return 3
        elif no_4>0 and no_2>0 and no_1==0:
            return 4
        else:
            print("unrecognised graph type, try increasing P")
            print("unrecognised graph: a=", a,"b= ", b)
            return 6
        
    def all_graph_type(self, n):
        '''compute graph type for all combos of coeffs less than n'''
        
        Z = np.empty((0,4), int)
        for a in range(-n,n):
            for b in range(-n,n):
                if (4*a**3 + 27*b**2) != 0:
                    Z = np.append(Z, [[a, b, no_int_roots(a,b), \
                            E_2.graph_type(self,a,b)]], axis=0) 
        return Z
    
    def all_graph_type_fast(self, n):
        '''compute graph type for all combos of coeffs less than n'''
        
        Z = np.empty((0,4), int)
        for a in range(-n,n):
            for b in range(-n,n):
                if (4*a**3 + 27*b**2) != 0:
                    Z = np.append(Z, [[a, b, no_int_roots(a,b), \
                            E_2.graph_type_fast(self,a,b)]], axis=0) 
        return Z
    
    def gen_graphtype_data(self, n):
        '''sort all graph type data into graph types for plotting of different colours'''
        
        Z = E_2.all_graph_type(self,n)
        types = [[],[],[],[],[],[]]
        for x in Z:
            if x[3] == 4 and x[2] == 0:
                types[4].append(list(x))
            else:
                types[x[3]-1].append(list(x))
        return types
    
    def gen_graphtype_data_fast(self, n):
        '''sort all graph type data into graph types for plotting of different colours'''
        
        Z = E_2.all_graph_type_fast(self,n)
        types = [[],[],[],[],[],[]]
        for x in Z:
            if x[3] == 4 and x[2] == 0:
                types[4].append(list(x))
            else:
                types[x[3]-1].append(list(x))
        return types
    
    def gen_graphtype_graph(self, n):
        '''plot the sorted graph type data for varied coefs'''
        
        Z = E_2.gen_graphtype_data(self, n)

        type1 = np.asarray(Z[0])
     #   type2 = np.asarray(Z[1])
        type3 = np.asarray(Z[2])
        type4 = np.asarray(Z[3])
#        type5 = np.asarray(Z[4])
      #  type6 = np.asarray(Z[5])
        
        plt.plot(type1[:,0], type1[:,1], 'y.')
      #  plt.plot(type2[:,0], type2[:,1], 'g.')
        plt.plot(type3[:,0], type3[:,1], 'r.')
        plt.plot(type4[:,0], type4[:,1], 'b.')
     #   plt.plot(type5[:,0], type5[:,1], 'm.')
    #    plt.plot(type6[:,0], type6[:,1], 'k.')
        plt.grid()
        plt.xlabel('a')
        plt.ylabel('b')
        plt.show()
    
    def gen_graphtype_graph_fast(self, n):
        '''plot the sorted graph type data for varied coefs'''
        
        Z = E_2.gen_graphtype_data_fast(self, n)

        type1 = np.asarray(Z[0])
     #   type2 = np.asarray(Z[1])
        type3 = np.asarray(Z[2])
        type4 = np.asarray(Z[3])
#        type5 = np.asarray(Z[4])
      #  type6 = np.asarray(Z[5])
        
        plt.plot(type1[:,0], type1[:,1], 'y.')
      #  plt.plot(type2[:,0], type2[:,1], 'g.')
        plt.plot(type3[:,0], type3[:,1], 'r.')
        plt.plot(type4[:,0], type4[:,1], 'b.')
     #   plt.plot(type5[:,0], type5[:,1], 'm.')
    #    plt.plot(type6[:,0], type6[:,1], 'k.')
        plt.grid()
        plt.xlabel('a')
        plt.ylabel('b')
        plt.show()
        
    def ratios_data(self, n):
        '''generate ratios between group sizes for ECs of graph type 1'''
        Z = np.empty((0,8), int)
        for a in range(-n,n):
            for b in range(-n,n):
                if (4*a**3 + 27*b**2) != 0:
                    no_elts = E_2.no_elts_ls_fast(self, a, b)
                    final = no_elts[-1]
                    if final[1]>0 and final[2]>0 and final[3]>0:
                        Z = np.append(Z, [[a, b, final[1], \
        final[2], final[3], final[1]/final[2], final[3]/final[1],\
            final[3]/final[2]]], axis=0)                      
        return Z
    
    def gen_cm_graph_ratios(self, n):
        '''plot cm for ratios between group sizes for ECs of graph type 1'''
        Z = E_2.ratios_data(self, n)
        #print(Z[0:20])
        x = np.reshape(Z[:, [0]], len(Z[:, [0]]))
        y = np.reshape(Z[:, [1]], len(Z[:, [1]]))
        z1 = np.reshape(Z[:, [5]], len(Z[:, 5]))
        z2 = np.reshape(Z[:, [6]], len(Z[:, 6]))
        z3 = np.reshape(Z[:, [7]], len(Z[:, 7]))
        plt.scatter(x, y, c=z1, cmap='coolwarm',vmin=min(z1), vmax=max(z1))
        plt.colorbar()
        plt.show()
        plt.scatter(x, y, c=z2, cmap='coolwarm',vmin=min(z2), vmax=max(z2))
        plt.colorbar()
        plt.show()
        plt.scatter(x, y, c=z3, cmap='coolwarm',vmin=min(z3), vmax=max(z3))
        plt.colorbar()
        plt.show()


gen_root_graph(20)
plt.show()
ec0 = E_2(2500)
ec0.plot_graph_no_elts(1,3)
ec0.plot_graph_props(1,5)
ec0.plot_graph_ratios(-3,1)
ec0.gen_graphtype_graph_fast(10)







