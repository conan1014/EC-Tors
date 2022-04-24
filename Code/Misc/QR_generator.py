import numpy as np
from sage.all import *
import time


N=10000

t0 = time.time()

def primes(P):
    '''generate all primes < P'''  
    return [i for i in range(2,P) if isPrime(i)]


def isPrime(a):
    '''fast prime checker'''
    return all(a % i for i in range(2, int(a)))


P = primes(N)
Z = np.empty((1), int)[1:]


for i in range(1,N+1):
    if i in P:
        Z = np.append(Z,[i],axis=0)
    else:
        Z = np.append(Z,[0],axis=0)
    
final = np.zeros((N,N), int)

for i in range(1,N+1):
    if Z[i-1]!=0:
        qrs = np.asarray(quadratic_residues(i))
        for x in qrs:
            final[i-1][x-1] = 1
    

print(final)
print(final[7-1][3-1])
t1 = time.time()

print(t1-t0)

#np.savetxt("foo.csv", final, delimiter=",")


