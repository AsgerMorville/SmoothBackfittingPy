# distutils: language = c++
# cython: language_level=3

cimport sbf
import numpy as np
import cython

def SBF_CPP(double[:] Y, double[:,:] X, double[:] output, int n, int d):
    #output_mat = np.zeros((dim1,dim2))
    #output_mat = output_mat.astype(float)
    #narr_view = cython.declare(cython.float[:, :, :], output_mat)
    # Directly pass the pointer to the data of the NumPy array
    #sumfunc.sumfunction1(&arr1[0,0],&arr2[0,0], &narr_view[0,0], dim1)
    sbf.sbfWrapper(&Y[0],&X[0,0], &output[0], n, d)
    return

def PL_SBF_CPP(double[:] Y, double[:,:] X, double[:,:] Z,  double[:] output, int n, int d_x, int d_z):
    sbf.plSBFWrapper(&Y[0],&X[0,0], &Z[0,0], &output[0], n, d_x, d_z)
    return