# Declare the class with cdef


cdef extern from "smooth_backfitting_core.h":
    void sbfWrapper(double* yPtr, double* xPtr, double* outputPtr, int n, int d)

cdef extern from "partially_linear_SBF.h":
    void plSBFWrapper(double* yPtr, double* xPtr, double* zPtr, double* outputPtr, int n, int d_x, int d_z)
