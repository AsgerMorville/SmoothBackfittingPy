# distutils: language = c++
# cython: language_level=3

cimport sbf
import numpy as np
import cython

def SBF_CPP(double[:] Y, double[:,:] X, double[:] output):
    n = X.shape[0]
    d = X.shape[1]
    #output_mat = np.zeros((dim1,dim2))
    #output_mat = output_mat.astype(float)
    #narr_view = cython.declare(cython.float[:, :, :], output_mat)
    # Directly pass the pointer to the data of the NumPy array
    #sumfunc.sumfunction1(&arr1[0,0],&arr2[0,0], &narr_view[0,0], dim1)
    sbf.sbfWrapper(&Y[0],&X[0,0], &output[0], n, d)
    return