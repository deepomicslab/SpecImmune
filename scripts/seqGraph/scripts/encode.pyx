import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)

cdef str _num_transfer(str seq):
    seq = seq.replace("A", "0").replace("C", "1").replace("G", "2").replace("T", "3")
    seq = ''.join(list(filter(str.isdigit, seq)))

    return seq    


cdef list _num_transfer_loc(str num_seq, int K):
    cdef list loc
    loc = []
    for i in range(0, len(num_seq)-K+1):
        loc.append(int(num_seq[i:i+K], 4))
    
    return loc

cdef np.ndarray[np.float64_t, ndim=1] _loc_transfer_matrix(list loc_list, int dis, int K):
    cdef np.ndarray[np.float64_t, ndim=2] matrix
    cdef np.ndarray[np.float64_t, ndim=1] new_matrix
    matrix = np.zeros((4**K, 4**K))
    for i in range(0, len(loc_list)-K-dis):
        matrix[loc_list[i]][loc_list[i+K+dis]] += 1

    new_matrix = matrix.flatten()
    
    return new_matrix

cdef np.ndarray[np.float64_t, ndim=1] _matrix_encoding(str seq, int K):
    cdef int length
    cdef np.ndarray[np.float64_t, ndim=1] feature
    cdef str num_seq
    cdef list loc
    seq = seq.upper()
    length = len(seq)
    num_seq = _num_transfer(seq)
    loc = _num_transfer_loc(num_seq, K)
    feature = np.hstack((
        _loc_transfer_matrix(loc, 0, K), 
        _loc_transfer_matrix(loc, 1, K), 
        _loc_transfer_matrix(loc, 2, K)))
    
    return feature/(length*1.0) * 100

def matrix_encoding(seq, K):

    return _matrix_encoding(seq, K)

