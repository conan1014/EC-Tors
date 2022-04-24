
from matplotlib import pyplot as plt
import numpy as np
from sage.all import *

'''small funcs'''

def sign(a):
    if a>=0:
        return '+'
    else:
        return ''

class ECWNF:
    """ECWNF: class of ECs in general WNF (y^2 = x^3 + ax^2 + bx + c)"""
    
    def __init__(self,*coeffs):
      self._c = list(coeffs)
      self.rps = []
      self.negrps = []
      self.linjoin = []
      self.vert = []
    
    def plot(self, n):
        """Plots the EC along with all stored rps, lines between rps and intermediary
        calculation points
        """
        
        '''
        # Alternative bounds of graph for greater scope
        listot = [item for sublist in self.rps for item in sublist]
        if len(listot)==0:
            x = np.arange(-10, 10, 0.1)
            y = np.arange(-10, 10, 0.1)
        else:
            coordmax = max(listot)
            delta = coordmax/500
            x = np.arange(-coordmax-1, coordmax+1, delta)
            y = np.arange(-coordmax-1, coordmax+1, delta)
        '''
        
        x = np.linspace(-n, n, 100)
        y = np.linspace(-n, n, 100)
        a, b, c = self._c[0], self._c[1], self._c[2]
        
        X, Y = np.meshgrid(x, y)
        Z = Y**2 - X**3 - self._c[0]*X**2 - self._c[1]*X - self._c[2]
        fig, ax = plt.subplots()
        CS = ax.contour(X, Y, Z, [0])
        plt.grid()
        plt.title("y^2 = x^3 {0} {1}x {2} {3}".format(sign(b), b, sign(c), c))
        for a in self.negrps:
            plt.plot(a[0], a[1], 'rx')
        for a in self.rps:
            plt.plot(a[0], a[1], 'bx')
        for a in self.linjoin:
            plt.plot(x, a[0]*x+a[1], '-r')
        for a in self.negrps:
            plt.vlines(a[0], a[1], -a[1],'r')
        for a in self.vert:
            plt.axvline(x=a, color='r')
        plt.axvline(x=0, color='r')
        
    def appendrps(self,*coords):
        """Add rp to EC"""
    
        self.rps.append(list(coords))
        
    def addrps(self, *rpsno):
        """Add together two rps using group addition, adding intemidiary calculation 
        points and lines to the EC
        """
    
        index = list(rpsno)

        x1, y1 = self.rps[abs(index[0])-1][0], np.sign(index[0])*self.rps[abs(index[0])-1][1]
        x2, y2 = self.rps[abs(index[1])-1][0], np.sign(index[1])*self.rps[abs(index[1])-1][1]
        

        a, b, c = self._c[0], self._c[1], self._c[2]
        if abs(index[0]) == abs(index[1]):
     
            if y1 == 0:
                self.vert.append(x1)
            else:
                print("duplication")
                x3 = (x1**4 - 2*b*x1**2 - 8*c*x1 - 4*a*c)/(4*x1**3 + 4*a*x1**2 + 4*b*x1 + 4*c)
                y3 = (x1**6 + 2*a*x1**5 + 5*b*x1**4 + 20*c*x1**3 + (20*a*c - 5*b**2)*x1**2 +\
                      (8*a**2*c - 2*a*b**2 - 4*b*c)*x1 + (4*a*b*c - b**3 -8*c**2))\
                    /(8*y1**3)
                m = (3*x1**2 + 2*a*x1 + b)/(2*y1)
                v = -m*x1 + y1
                self.rps.append([x3,y3])
                self.negrps.append([x3,-y3])
                self.linjoin.append([m,v])
            
        else:
            print("normal addition")
            m = (y2-y1)/(x2-x1)
            v = y1 - m*x1
            x3 = m**2 - a - x1 - x2
            y3 = m*x3 + v
            self.rps.append([x3,-y3])
            self.negrps.append([x3,y3])
            self.linjoin.append([m,v])
            
    def ord2(self):
        """ Generate all complex, and plot all rational points of order 2
        """
        
        coeff = [1, self._c[0], self._c[1], self._c[2]]
        roots = np.roots(coeff)
        real_roots = [x.real for x in roots if x.imag==0]
        for x in real_roots:
            self.rps.append([x,0])
            self.negrps.append([x,0])
        return(roots)
    
    def ord3(self):
        """ Generate all complex, and plot all rational points of order 3
        THIS METHOD IS UNFINISHED
        """
        
        a, b, c = self._c[0], self._c[1], self._c[2]
        coeff = [3, 4*a, 6*b, 12*c, (4*a*c-b**2)]
        xroots_temp = np.roots(coeff)
        roots = []
        for x in xroots_temp:
            y = np.sqrt(x**3 + a*x**2 + b*x + c)
            roots.append([x,y])
            roots.append([x,-y])
        real_roots = [[x[0].real,x[1].real] for x in roots if x[0].imag==0 \
                      and x[1].imag==0]
        for x in real_roots:
            self.rps.append(x)
            self.negrps.append([x[0],-x[1]])
        return(roots)
    
    def finite_ord_poss(self):
        """ Generate and plot all possible finite order rational points
        THIS METHOD IS UNFINISHED
        """
        
        a, b, c = self._c[0], self._c[1], self._c[2]
        D = abs(-4*a**3+a**2*b**2 + 18*a*b*c - 4*b**3 - 27*c**2)
        y_poss = [int(np.sqrt(i)) for i in range(1, D + 1) if \
                    D % i == 0 and np.sqrt(i).is_integer()]
        final = []
        for y in y_poss:

            coeff = [1, a, b, c-y**2]
            x1 = np.roots(coeff)
  
            x2 = [x.real for x in x1 if x.imag==0]

            x3 = [round(x,3) for x in x2 if round(x,3).is_integer()]
      
            for x in x3:
                final = final + [[x, y]]
                final = final + [[x, -y]]
                
        coeff0 = [1, self._c[0], self._c[1], self._c[2]]
        roots0 = np.roots(coeff0)
        real_roots0 = [x.real for x in roots0 if x.imag==0]
        int_roots0 = [round(x,3) for x in real_roots0 if round(x,3).is_integer()]
        
        for x in int_roots0:
            final = [[x,0]] + final
        
        
        
        for x in final:
            self.rps.append(x)
            self.negrps.append([x[0],-x[1]])
        
        return(final)
    





ec1 = ECWNF(0,3,4)
ec1.appendrps(0,2)
ec1.appendrps(0,-2)
ec1.addrps(1,2)
ec1.plot(5)

