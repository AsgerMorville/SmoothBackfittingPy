# distutils: language = c++
# cython: language_level=3

cimport sbf
import numpy as np
import cython

def SBF_CPP(double[:] Y, double[:,:] X, double[:] output, int n, int d):
    sbf.sbfWrapper(&Y[0],&X[0,0], &output[0], n, d)
    return

def PL_SBF_CPP(double[:] Y, double[:,:] X, double[:,:] Z,  double[:] output, int n, int d_x, int d_z):
    sbf.plSBFWrapper(&Y[0],&X[0,0], &Z[0,0], &output[0], n, d_x, d_z)
    return

def SBF(Y, X):
    # Make sure theyre stored contiguously in memory
    Y = np.ascontiguousarray(Y)
    X = np.ascontiguousarray(X)
    n, d = X.shape
    output_vec = np.zeros(n)
    SBF_CPP(Y, X, output_vec, n, d)
    return output_vec

def PL_SBF(Y, X, Z):
    # Make sure theyre stored contiguously in memory
    Y = np.ascontiguousarray(Y)
    X = np.ascontiguousarray(X)
    Z = np.ascontiguousarray(Z)
    n, d_x = X.shape
    d_z = Z.shape[1]
    output_vec = np.zeros(n)
    PL_SBF_CPP(Y, X, Z, output_vec, n, d_x, d_z)
    return output_vec

