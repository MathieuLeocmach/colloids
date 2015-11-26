import numpy as np
cimport numpy as np
from libc.math cimport floor, sqrt

DTYPE = np.int
ctypedef np.int_t DTYPE_t

def generate_kvec3D(int nk, np.ndarray[DTYPE_t, ndim=2] kvectors):
    """Generate the indices of the 3D wavevectors in a cubic box"""
    cdef int maxNvec = len(kvectors)
    cdef int i,j,k
    cdef int nvec = 0
    for i in range(nk):
        if nvec >= maxNvec: break
        #for j in range(int(sqrt(nk*nk-i*i))):
        j=0
        while j**2 < nk**2 - i**2:
            if nvec >= maxNvec: break
            k = int(sqrt(nk*nk-i*i-j*j))
            #if k != int(k): continue
            kvectors[nvec,0] = i
            kvectors[nvec,1] = j
            kvectors[nvec,2] = k
            kvectors[nvec+1,0] = j
            kvectors[nvec+1,1] = k
            kvectors[nvec+1,2] = i
            kvectors[nvec+2,0] = k
            kvectors[nvec+2,1] = i
            kvectors[nvec+2,2] = j
            nvec += 3
            j += 1
    return kvectors[:maxNvec]
