import os
import math 
import numpy as np
from numpy import linalg as la

numSam = 10
numEval = 1000
numSiteX = 6
numSiteY = 2
mz = 0.5

numSite = numSiteX*numSiteY

# Contruct the basis of the Hilbert space.
max = int(math.pow(2, numSite-1)) # maxmium dimension of the Hilbert space for the sub-space in terms of pure spin configuration
n = 0  # counter for states in a specific magnetic sector
spinBasis = np.zeros(0, dtype=int)
for i in range(max):
    if bin(i).count('1') == int(mz+(numSite-1)/2.0):
        # print(bin(i))
        spinBasis = np.append(spinBasis, i)
        n += 1
holeBasis = np.zeros(numSite, dtype=int)
for i in range(numSite):
    holeBasis[i] = i

subDim = n  # dimension of the sub-space with a hole fixed
dim = numSite*n  # dimension of the Hilbert space

def LoadEigVal(f):
    rawData = np.fromfile(f, dtype=np.float)
    eigVal = np.reshape(rawData, (numSam, numEval))
    return eigVal

# The returned eigVec[i][j] gives the ith sample's jth eigvector with dimesion dim.
def LoadEigVec(f):
    rawData = np.fromfile(f, dtype=np.float)
    temp = np.reshape(rawData, (numSam, numEval, dim, 2))
    eigVec = np.array(temp[..., 0], dtype=complex)
    eigVec.imag = temp[..., 1]
    # print(eigVec)
    return eigVec

def FindState(s, lower, upper):
    if lower == upper:
        return upper
    else:
        mid = (lower+upper)//2
        if s > spinBasis[mid]:
            return FindState(s, mid+1, upper)
        else:
            return FindState(s, lower, mid)

def Flip(s, i, j):
    f = int(np.power(2, i))+int(np.power(2, j))
    return s^f

# Translation() is only valid for periodic boundary condition. Note that there is also a fermion sign arised here if the hole is not located in N-1.

def TranslationX(v):
    w = np.zeros(dim, dtype=complex)
    for l in range(dim):
        h = holeBasis[l // subDim] # h = hy*numSiteX+hx 
        hx = h % numSiteX
        hy = h // numSiteX
        hh = hy*numSiteX+((hx+1) % numSiteX)

        s = spinBasis[l % subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'0'+b[h:] # Reconstruct the half-filled configuration.
        lst = [None]*numSite
        for j in range(numSiteY):
            for i in range(numSiteX):
                k = j*numSiteX+i
                kk = j*numSiteX+((i+1) % numSiteX)
                lst[kk] = config[k]
        newConfig = ''.join(lst)
        sign = np.power(-1., hh-h-(numSiteX-1)*numSiteY)  # fermionic sign arised  
        bb = newConfig[:hh]+newConfig[hh+1:]
        bb = bb[::-1]
        ss = int(bb, 2)
        n = FindState(ss, 0, subDim-1)
        w[hh*subDim+n] = sign*v[l]
    return w

def TranslationY(v):
    w = np.zeros(dim, dtype=complex)
    for l in range(dim):
        h = holeBasis[l // subDim] # h = hy*numSiteX+hx 
        hx = h % numSiteX
        hy = h // numSiteX
        hh = ((hy+1) % numSiteY)*numSiteX+hx

        s = spinBasis[l % subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'0'+b[h:] # Reconstruct the half-filled configuration.
        lst = [None]*numSite
        for j in range(numSiteY):
            for i in range(numSiteX):
                k = j*numSiteX+i
                kk = ((j+1) % numSiteY)*numSiteX+i
                lst[kk] = config[k]
        newConfig = ''.join(lst)
        # print(h, hh, config, newConfig)
        sign = np.power(-1.0, hh-h)
        bb = newConfig[:hh]+newConfig[hh+1:]
        bb = bb[::-1]
        ss = int(bb, 2)
        n = FindState(ss, 0, subDim-1)
        w[hh*subDim+n] = sign*v[l]
    return w

def Rotation( v):
    # Generate the permutation for couterclockwise rotation pi/2.
    per = [None]*numSite
    for j in range(numSiteY):
        for i in range(numSiteX):
            k = j*numSiteX+i
            kk = i*numSiteY+(numSiteX-1-j)
            per[k] = kk
    # print(per)

    w = np.zeros(dim, dtype=complex)
    for l in range(dim):
        h = holeBasis[l // subDim] # h = hy*numSiteX+hx 
        hh = per[h]

        s = spinBasis[l % subDim]
        b = bin(s)[2:]
        b = b.rjust(numSite-1, '0')
        b = b[::-1]
        config = b[:h]+'0'+b[h:] # Reconstruct the half-filled configuration.
        lst = [None]*numSite
        for j in range(numSiteY):
            for i in range(numSiteX):
                k = j*numSiteX+i
                kk = i*numSiteY+(numSiteX-1-j)
                lst[per[k]] = config[k]
        newConfig = ''.join(lst)
        sign = np.power(-1.0, hh-h)  # fermionic sign arised  
        bb = newConfig[:hh]+newConfig[hh+1:]
        bb = bb[::-1]
        ss = int(bb, 2)
        n = FindState(ss, 0, subDim-1)
        w[hh*subDim+n] = sign*v[l]
    return w

def MarshallSign(i):
    s = spinBasis[i % subDim]
    h = holeBasis[i // subDim]
    b = bin(s)[2:]
    b = '%0*.0f' % (numSite, int(b))
    b = b[::-1]
    n = 0
    for j in range(numSite):
        if 0 == j % 2:
            if j == h:
                n += 1
            elif j < h and '0' == b[j]:
                n += 1
            elif j > h and '0' == b[j-1]:
                n += 1
    return np.power(-1, n)
