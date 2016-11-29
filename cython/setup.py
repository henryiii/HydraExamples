from Cython.Build import cythonize
from distutils.core import setup

setup(
        name='Hello world app',
        ext_modules = cythonize('Hydra.pyx'),
        )
