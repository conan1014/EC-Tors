from ECWNF_class.py import *
import math
from scipy.integrate import quad
import matplotlib.pyplot as plt

class s_t_comparison:
    '''class to verify the Sato-Tate conjecture'''
    
    def __init__(self, P):
        self.P = P
        self.biglist = s_t_comparison.thetap(self)
    
    def is_prime(m):
        
        f = 2
        while True and f*f <= m:
            if m%f == 0:
                return False
            f = f+1
        return True
    
    def Mp(self):  
        
        '''generate Mp (size of E(F_p)) for all primes less than P'''
        
        lst1 = [i for i in range(self.P) if s_t_comparison.is_prime(i)]
        lst2 = [[None, None]]
        for x in lst1:
            ec1 = EC_modp(0,1,1,x)
            Mp = len(ec1.genpoints())
            lst2 = np.append(lst2, [[x, Mp]], axis=0)
        return lst2

    def thetap(self):
        '''generate theta_p for all primes less than P'''
        
        lst1 = [i for i in range(self.P)[1:] if s_t_comparison.is_prime(i)]
        lst2 = [[None,None,None]]
        for x in lst1:
            ec1 = EC_modp(0,1,1,x)
            Mp = len(ec1.genpoints())
            thetap = math.acos((Mp-x-1)/(2*np.sqrt(x)))
            lst2 = np.append(lst2, [[x, Mp, round(thetap,3)]], axis=0)
        return lst2[1:]

    def plot_Mp(self):
        '''Plot Mp by p'''
        
        lst = s_t_comparison.Mp(self)
        X = lst[:,0]
        Y = lst[:,1]
        plt.plot(X,Y)
    
    def plot_thetap(self):
        '''Plot theta_p by p'''
        
        lst = s_t_comparison.thetap(self)
        X = lst[:,0]
        Y = lst[:,2]
        
        plt.plot(X,Y)
    
    def vals_bet(self, a,b):
        '''returns number of points in Mp between a and b/no of primes less than P'''
        
        lst = self.biglist
        lst_red = [x for x in lst[1:] if a<= x[2] <=b]
        return len(lst_red)/len(lst)
    
    def integrand(t):
        return (math.sin(t))**2
    
    def integral(a,b):
        '''returns integral calculation'''
        
        result = quad(s_t_comparison.integrand,a,b)
        return (2/math.pi)*result[0]
    
    def st_comparison(self, a,b):
        '''returns comparison'''
        
        return round(s_t_comparison.vals_bet(self, a,b),3), round(s_t_comparison.integral(a,b),3)
  


st1 = s_t_comparison(200)
st1.plot_thetap()
print(st1.st_comparison(0, math.pi/2))
print(st1.st_comparison(math.pi/3, math.pi/2))
print(st1.st_comparison(2*math.pi/5, 3*math.pi/5))
























