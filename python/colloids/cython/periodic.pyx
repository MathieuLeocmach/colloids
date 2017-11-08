import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True) # turn off division by zero checking for entire function
def periodify(
    np.ndarray[np.float64_t, ndim=2] u, 
    np.ndarray[np.float64_t, ndim=2] v, 
    periods=None
):
    assert u.shape[0] == v.shape[0]
    assert u.shape[1] == v.shape[1]
    #cdef np.ndarray[np.float64_t, ndim=2] diff = v - u
    if periods is None:
        return v - u
    assert len(periods) == u.shape[1]
    cdef np.ndarray[np.float64_t, ndim=2] diff = np.zeros_like(u)
    #ensures the largest coordinate is smaller than half a period
    cdef np.ndarray[np.float64_t, ndim=1] pers = np.array(periods)
    cdef np.ndarray[np.float64_t, ndim=1] half = 0.5*np.array(periods)
    cdef int i,j
    cdef float d,h,p
    for i in range(diff.shape[0]):
        for j in range(diff.shape[1]):
            d = v[i,j] - u[i,j]
            h = half[j]
            p = pers[j]
            while d > h:
                d -= p
            while d < -h:
                d += p
            diff[i,j] = d
    return diff
