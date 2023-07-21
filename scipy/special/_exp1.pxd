cdef extern from "_exp1.cpp" nogil:
    double _exp1(double)


cdef inline double exp1(double x) nogil:
    return _exp1(x)
