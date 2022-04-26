import numpy as np


gens_gl23 = [[2,1,1,1], [1,1,1,2], [1,2,0,1], [2,0,0,1]]

gens_d12 = [[2,0,0,1], [2,1,0,1], [1,2,0,1], [1,0,0,2]]

gens_sd16 = [[0,1,1,1], [1,1,0,2]]

gens_s3 = [[1,1,0,1], [1,2,0,2]]

gens_c2 = [[1,0,0,1], [1,0,0,2]]

gens_d4 = [[0,1,1,0], [1,0,0,2]]

gens_c2c2 = [[1,0,0,2], [2,0,0,1]]




def is_in(gp, temp):
    for el in gp:
        if el == temp:
            return True
    else:
        return False
    
def matmul(a,b,p):
    '''matrix product modulo p'''
    A = np.array([[a[0], a[1]],[a[2],a[3]]])
    B = np.array([[b[0], b[1]],[b[2],b[3]]])
    C = np.matmul(A,B)
    mod = np.array([[p,p],[p,p]])
    C = np.mod(C, mod)
    return [C[0][0], C[0][1], C[1][0], C[1][1]]

def inv(a, gp, p):
    '''matrix inverse'''
    for b in gp:
        if matmul(a,b,p)==[1,0,0,1]:
            return b

def group(gens, size, p):
    '''generate the whole group from generators'''
    gp = gens
    for a in gp:
        for b in gp:
            c = matmul(a,b,p)
            if not is_in(gp,c):
                gp.append(c)
    if len(gp)==size:
        return gp
    else:
        return group(gp,size,p)
    
def ccl(x, gp, p):
    '''generate the conjugacy class of one element'''
    ccl = []
    for el in gp:
        el_inv = inv(el, gp, p)
        temp = matmul(matmul(el, x, p), el_inv, p)
        if not is_in(ccl, temp):
            ccl.append(temp) 
    return ccl

def ccls(gp,p):
    '''generate all conjugacy classes'''
    ccls = []
    ccls.append(ccl(gp[0], gp, p))
    for x in gp:
        already_there = False
        for cl in ccls:
            if is_in(cl, x):
                already_there = True
        if not already_there:
            ccls.append(ccl(x, gp, p))

    return ccls




gp = group(gens_gl23,48,3)
for cl in ccls(gp, 3):
    print(cl)
    print()

        

    
























