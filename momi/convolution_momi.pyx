cimport numpy as np
import numpy as np
import cython


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def convolve_momi(np.ndarray[np.double_t, ndim=2] A,
                  np.ndarray[np.double_t, ndim=2] B):
    # first dimension is for snps (shared by A and B)
    assert A.shape[0] == B.shape[0]
    # second dimension is for derived allele counts
    
    cdef np.ndarray[np.double_t, ndim=2] C = np.zeros((A.shape[0],
                                                       A.shape[1] + B.shape[1] - 1))

    cdef unsigned int i,l,m,I,L,M
    I,L = A.shape[0],A.shape[1]
    M = B.shape[1]

    for i in range(I):
        for l in range(L):
            for m in range(M):
                C[i,l+m] += A[i,l] * B[i,m]
    return C

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def convolve_chen(np.ndarray[np.double_t, ndim=3] A,
                  np.ndarray[np.double_t, ndim=3] B):
    # first dimension is for snps (shared by A and B)
    assert A.shape[0] == B.shape[0]
    # second, third dimensions are for ancestral/derived allele counts
    assert A.shape[1] == A.shape[2] and B.shape[1] == B.shape[2]
    
    cdef np.ndarray[np.double_t, ndim=3] C = np.zeros((A.shape[0],
                                                       A.shape[1] + B.shape[1] - 1,
						       A.shape[2] + B.shape[2] - 1))

    cdef unsigned int i,l,ll,m,mm,I,L,M
    I,L = A.shape[0],A.shape[1]
    M = B.shape[1]

    for i in range(I):
        for l in range(L):
            for ll in range(L-l):
                for m in range(M):
                    for mm in range(M-m):
                        C[i,l+m,ll+mm] += A[i,l,ll] * B[i,m,mm]
    return C
