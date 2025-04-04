# cython: language_level=3, boundscheck=False, wraparound=False

import numpy as np
cimport numpy as cnp
ctypedef cnp.float64_t DTYPE_T #numpyの詳細な型指定

class Example2:
    
    def __init__(self, cnp.ndarray[DTYPE_T, ndim=1] A, cnp.ndarray[DTYPE_T, ndim=1] B):

        self.A = A
        self.B = B

    def sum_of_array(self):

        # cdef cnp.ndarray[double, ndim=1] A = self.A
        # cdef cnp.ndarray[double, ndim=1] B = self.B
        cdef cnp.ndarray[DTYPE_T, ndim=1] result

        result = self.A + self.B

        return result