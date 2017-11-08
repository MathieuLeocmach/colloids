from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Colloid toolkit',
  ext_modules = cythonize("periodic.pyx"),
)
