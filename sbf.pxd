# Declare the class with cdef


cdef extern from "smooth_backfitting_core.h":
    void sbfWrapper(double* yPtr, double* xPtr, double* outputPtr, int n, int d)

