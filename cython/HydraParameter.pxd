from libcpp.string cimport string

cdef extern from "Parameter.h" namespace "hydra":
    cdef cppclass Parameter:
        Parameter() except +
        Parameter(string name, double value, double error, double downlim, double uplim ) except +
        Parameter(string name, double value, double error) except +
        Parameter(string name, double value) except +

