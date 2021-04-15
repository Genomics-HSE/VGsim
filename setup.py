from setuptools import setup

from Cython.Build import cythonize
from Cython.Compiler import Options
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy

Options.annotate = True

ext = Extension("_BirthDeath", ["_BirthDeath.pyx"],
                include_dirs = [numpy.get_include()],
                language='c++',)

setup(ext_modules=[ext],
      cmdclass = {'build_ext': build_ext})
