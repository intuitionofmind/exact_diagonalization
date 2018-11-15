import math 
import numpy as np
from numpy import linalg as la

numSite = 10
numEval = 1200
numSam = 5
step = 0.1
mz = 0.5  # total magnetization, a specific sector
flux = 0.0
cut = 5

inter = 10
fs = 16
ls = 12

max = int(math.pow(2, numSite-1)) # maxmium dimension of the Hilbert space for the sub-space in terms of pure spin configuration
n = 0  # counter for states in a specific magnetic sector
spinBasis = np.zeros(0, dtype=int)
for i in range(max):
    if bin(i).count('1') == int(mz+(numSite-1)/2.0):
        spinBasis = np.append(spinBasis, i)
        n += 1

holeBasis = np.zeros(numSite, dtype=int)
for i in range(numSite):
    holeBasis[i] = i

subDim = n  # dimension of the sub-space with a hole fixed
dim = numSite*n  # dimension of the Hilbert space

print(subDim, dim)

fe = '/Users/wayne/Downloads/data/eigenvalues%s_size1D%s_J0.1_0.1_%s_%s_sigma%s.dat'
fv = '/Users/wayne/Downloads/data/eigenvectors%s_size1D%s_J0.1_0.1_%s_%s_sigma%s.dat'

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

# Hamiltonian() compute t-J and \sigma t-J spin chain with periodic boundary condition.

def Hamiltonian(v, J, sigma):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        s = spinBasis[i % subDim] # i = h*sumDim+(i % subDim)
        h = holeBasis[i // subDim]
        b = bin(s)[2:]
        b.rjust(numSite, '0')
        b = b[::-1]
        for j in range(numSite):
            sign = 1.0
            jj = (j+1) % numSite
            if j != h and jj != h:
                if j < h and jj < h and b[j] != b[jj]:
                    w[i] -= 0.5*J*v[i]
                    ss = Flip(s, j, jj)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*J*v[i]
                elif j > h and jj > h and b[j-1] != b[jj-1]:
                    w[i] -= 0.5*J*v[i]
                    ss = Flip(s, j-1, jj-1)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*J*v[i]
                elif jj == 0 and j == (numSite-1) and b[j-1] != b[jj]: # Only for the case jj == 0 and j == numSite-1 of periodic boundary condition
                    w[i] -= 0.5*J*v[i]
                    ss = Flip(s, j-1, jj)
                    l = FindState(ss, 0, subDim-1)
                    w[h*subDim+l] += 0.5*J*v[i]
            elif jj == h:
                hh = (h-1+numSite) % numSite
                if 0 == h: # This is a special case in which the hole hopping across the boundary and the spin configuration will also be changed.
                    if sigma and '0' == b[j-1]:
                        sign = -1.0
                    bb = b[numSite-2]+b[:numSite-2]
                    bb = bb[::-1]
                    ss = int(bb, 2)
                    l = FindState(ss, 0, subDim-1)
                    w[hh*subDim+l] -= 1.0*np.power(-1, numSite)*sign*v[i]
                else:
                    if sigma and '0' == b[j]:
                        sign = -1.0
                    w[hh*subDim+(i % subDim)] -= 1.0*sign*v[i]
            elif j == h:
                hh = (h+1) % numSite
                if (numSite-1) == h:
                    if sigma and '0' == b[jj]:
                        sign = -1.0
                    bb = b[1:numSite-1]+b[0]
                    bb = bb[::-1]
                    ss = int(bb, 2)
                    l = FindState(ss, 0, subDim-1)
                    w[hh*subDim+l] -= 1.0*np.power(-1, numSite)*sign*v[i]
                else:
                    if sigma and '0' == b[jj-1]:
                        sign = -1.0
                    w[hh*subDim+(i % subDim)] -= 1.0*sign*v[i]
    return w

def PrintHam(J, sigma):
    H = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        v1 = np.zeros(dim)
        v1[i] = 1.0
        for j in range(dim):
            v2 = np.zeros(dim)
            v2[j] = 1.0
            vv2 = Hamiltonian(v2, J, sigma)
            H[i][j] = np.vdot(v1, vv2)
    print(np.real(H))

# Translation() is only valid for periodic boundary condition. Note that there is also a fermion sign arised here if the hole is not located in N-1.

def Translation(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        h = holeBasis[i // subDim]
        hh = (h+1) % numSite
        if h == numSite-1:
            w[hh*subDim+(i % subDim)] = v[i]
        else:
            s = spinBasis[i % subDim]
            b = bin(s)[2:]
            b = b.rjust(numSite, '0')
            b = b[::-1]
            bb = b[numSite-2]+b[:numSite-2]
            bb = bb[::-1]
            ss = int(bb, 2)
            l = FindState(ss, 0, subDim-1)
            w[hh*subDim+l] = np.power(-1.0, numSite)*v[i]
    return w

'''
def Translation(v):
    w = np.zeros(dim, dtype=complex)
    for i in range(dim):
        h = holeBasis[i // subDim]
        hh = (h+1) % numSite
        if h == numSite-1:
            w[hh*subDim+(i % subDim)] = v[i]
        else:
            s = spinBasis[i % subDim]
            b = bin(s)[2:]
            b = b.rjust(numSite, '0')
            b = b[::-1]
            bb = b[numSite-2]+b[:numSite-2]
            bb = bb[::-1]
            ss = int(bb, 2)
            l = FindState(ss, 0, subDim-1)
            w[hh*subDim+l] = np.power(-1.0, numSite)*v[i]
    return w
'''
 
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

# ---------------------Test begin--------------------

'''
PrintHam(1.2, True)
print(test_paras)
ene = LoadEigVal(test_paras)
wf = LoadEigVec(test_paras, 0)
wfT = LoadTransEigVec(test_paras, 0)

v = wf[slice]
print('Check energy: ', ene[slice][0], np.vdot(v, Hamiltonian(v, (slice+1)*0.1, False)))
vv = Translation(v)
for i in range(dim):
    print(vv[i]-wfT[slice][i])
'''

'''
# Test Translation() and diagonalization.

ene = LoadEigVal(test_paras)
wf0 = LoadEigVec(test_paras, 0)
wf1 = LoadEigVec(test_paras, 1)
temp = np.zeros((2, 2), dtype=complex)
if np.absolute(ene[slice][0]-ene[slice][1]) < 0.0000000001:
    wf0t = Translation(wf0[slice])
    wf1t = Translation(wf1[slice])
    temp[0][0] = np.vdot(wf0[slice], wf0t) # vdot for complex conjugate...
    temp[0][1] = np.vdot(wf0[slice], wf1t)
    temp[1][0] = np.vdot(wf1[slice], wf0t)
    temp[1][1] = np.vdot(wf1[slice], wf1t)
    w, v = la.eig(temp)
    print(slice, np.angle(w, deg=True))
    if w[0].imag > 0: # choose the eige
        choice = 0
    else:
        choice = 1
    wf = v[:, choice][0]*wf0[slice]+v[:, choice][1]*wf1[slice]  #construct eigenstates with specific momentum
    wff = Translation(wf)
    print('Test trsanlation:')
    for i in range(10):
        print(np.angle(wf[i]/wff[i], deg=True))
'''
# ----------------------Test end----------------------
