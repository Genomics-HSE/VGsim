@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose1(double[::1] w, double tw, double rn):
    cdef:
        Py_ssize_t i
        double total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
        for i in range(w.shape[0]):
            print_err(w[i], end=", ")
        print_err()
        print_err(tw)
        print_err(rn)
        sys.exit(1)
    return [ i, ( rn-(total-w[i]) )/w[i] ]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose2(long[::1] w, long tw, double rn):
    cdef:
        Py_ssize_t i
        double total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose1_skip(double[::1] w, double tw, double rn, Py_ssize_t skip):
    cdef:
        Py_ssize_t i
        double total

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
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose2_skip(long[::1] w, long tw, double rn, Py_ssize_t skip):
    cdef:
        Py_ssize_t i
        double total

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
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]
