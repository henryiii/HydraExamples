# distutils: language = c++

cimport HydraParameter as hp
from libcpp.string cimport string

def HydraHello():
    return "Hello from Cython"

cdef class Parameter:
    cdef hp.Parameter* thisptr
    def __cinit__(self, string name, double value):
        self.thisptr = new hp.Parameter(name, value)

    def __dealloc__(self):
        del self.thisptr
