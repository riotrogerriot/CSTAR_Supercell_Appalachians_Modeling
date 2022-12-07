from distutils.core import setup
from Cython.Build import cythonize
setup(name="srh",ext_modules=cythonize('srh_cy.pyx'))
