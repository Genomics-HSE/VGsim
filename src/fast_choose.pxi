ctypedef fused double_or_npy_int64:
    double
    npy_int64

cdef inline void print_error(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn, str func_name):
    cdef Py_ssize_t i
    print(f"{func_name} alert: 0-weight sampled")
    for i in range(w.shape[0]):
        print(w[i], end=", ")
    print()
    print(tw)
    print(rn)
    sys.exit(1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn):
    cdef:
        Py_ssize_t i
        double_or_npy_int64 total

    rn = tw * rn
    i = 0
    total = w[0]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_error(w, tw, rn, "fastChoose()")
    return [i, (rn - (total - w[i])) / w[i]]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose_skip(double_or_npy_int64[::1] w, double_or_npy_int64 tw, double rn, Py_ssize_t skip):
    cdef:
        Py_ssize_t i
        double_or_npy_int64 total

    rn = tw * rn
    i = 0
    if skip == 0:
        i += 1
    total = w[i]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        if i != skip:
            total += w[i]
    if w[i] == 0.0:
        print_error(w, tw, rn, "fastChoose_skip()")
    return [i, (rn - (total - w[i])) / w[i]]
