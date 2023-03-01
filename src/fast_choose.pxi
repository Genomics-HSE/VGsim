import numpy as np
def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

ctypedef fused double_or_npy_int64:
    double
    npy_int64


cdef inline void print_error(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn):
    cdef Py_ssize_t i
    for i in range(w.shape[0]):
        print_err(w[i], end=", ")
    print_err()
    print_err(tw)
    print_err(rn)
    sys.exit(1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn):
    cdef:
        Py_ssize_t i
        double_or_npy_int64 total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print("fastChoose() alert: 0-weight sampled")
        print_error(w, tw, rn)
    return [ i, ( rn-(total-w[i]) )/w[i] ]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose_vec(vector[double] w, double tw, double rn):
    cdef:
        Py_ssize_t i
        double total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.size() - 1:
        i += 1
        total += w[i]
    #if w[i] == 0.0:
       # print("fastChoose_vec() alert: 0-weight sampled")
        # print_error(w, tw, rn)
    return [ i, ( rn-(total-w[i]) )/w[i] ]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose_skip(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn, Py_ssize_t skip):
    cdef:
        Py_ssize_t i
        double_or_npy_int64 total

    rn = tw*rn
    if skip > 0:
        i = 0
        total = w[0]
    else:
        i = 1
        total = w[1]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        if i == skip:
            i += 1
        total += w[i]
    if w[i] == 0.0:
        print("fastChoose_skip() alert: 0-weight sampled")
        print_error(w, tw, rn)
    return [ i, ( rn-(total-w[i]) )/w[i] ]
