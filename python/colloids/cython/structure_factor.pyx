cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport floor, sqrt, sin, M_PI
cdef extern from "complex.h":
    double complex cexp(double complex)
    double creal(double complex);
    double cimag(double complex);

ctypedef np.int_t DTYPE_t

@cython.boundscheck(False)
cdef DTYPE_t generate_kvec3D(int nk, np.ndarray[DTYPE_t, ndim=2] kvectors):
    """Generate the indices of the 3D wavevectors in a cubic box"""
    cdef DTYPE_t maxNvec = kvectors.shape[0]
    assert maxNvec%3 == 0
    cdef int i,j,k
    cdef DTYPE_t nvec = 0
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
    return nvec

@cython.boundscheck(False)
cdef DTYPE_t generate_kvec2D(int nk, np.ndarray[DTYPE_t, ndim=2] kvectors):
    """Generate the indices of the 2D wavevectors in a square box"""
    cdef DTYPE_t maxNvec = kvectors.shape[0]
    cdef double M = 0.5*M_PI/maxNvec
    cdef int i,j,u
    cdef DTYPE_t nvec = 0
    
    for u in range(maxNvec/2):
        i = int(nk * sin(M * u))
        if nvec >0 and kvectors[nvec-1, 0] == i:
            continue
        j = int(sqrt(nk**2 - i**2))
        kvectors[nvec,0] = i
        kvectors[nvec,1] = j
        kvectors[nvec+1,0] = j
        kvectors[nvec+1,1] = i
        nvec += 2

    return nvec

def get_Sq(np.ndarray[np.float64_t, ndim=2] pos, int Nbins, double L=203.0, DTYPE_t maxNvec=30, np.ndarray[double complex, ndim=1]field=None):
    """Structure factor of the positions in a square or cubic box of size L."""
    #implementation without parallelism, 
    #within 10% of previous scipy.weave implementation (parallelism disabled)
    if pos.shape[1] ==2:
        generate_kvec = generate_kvec2D
    elif pos.shape[1] == 3:
        generate_kvec = generate_kvec3D
    else:
        raise NotImplementedError("generation of k-vectors in {} dimensions is not implemented".format(pos.shape[1]))
    if field is None:
        field = np.ones(pos.shape[0], np.complex128)
    assert field.shape[0] == pos.shape[0]
    cdef double complex prod, sum_rho, e
    cdef DTYPE_t i, j, k, nk, kv, nvec
    cdef double p
    cdef np.ndarray[DTYPE_t, ndim=2] kvectors= np.zeros((maxNvec, pos.shape[1]), dtype=int)
    cdef np.ndarray[double, ndim=1] Sq = np.zeros(Nbins)
    #cache the complex exponential caculations
    cdef np.ndarray[double complex, ndim=3] cache = np.ones([Nbins, pos.shape[0], pos.shape[1]], np.complex128)
    for i in range(pos.shape[0]):
        for j in range(pos.shape[1]):
            p = pos[i,j]
            e = cexp(2j*M_PI/L * p)
            for nk in range(1,Sq.shape[0]):
                cache[nk, i, j] = cache[nk-1, i, j] * e
    for nk in range(1,Sq.shape[0]):
        #generate the k-vectors
        nvec = generate_kvec(nk, kvectors)
        for k in range(nvec):
            sum_rho = 0j
            for i in range(pos.shape[0]):
                prod = field[i]
                for j in range(pos.shape[1]):
                    kv = kvectors[k,j]
                    #random access to cache
                    prod *= cache[kv, i, j]
                sum_rho += prod
            Sq[nk] += creal(sum_rho)**2 + cimag(sum_rho)**2
        Sq[nk] /= nvec
    return Sq / pos.shape[0]

    

