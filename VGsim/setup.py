from setuptools import setup

from Cython.Compiler import Options

import sys
import numpy

Options.annotate = True

PACKAGE_NAME = 'VGsim'

# enable openMP on linux, disable on MacOS
if sys.platform == "darwin":
    openmp_args = []
else:
    openmp_args = ['-fopenmp']


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(PACKAGE_NAME, parent_package, top_path)

    config.add_extension('_BirthDeath',
                         sources=['_BirthDeath.cxx'],
                         depends=['_BirthDeath.pyx',
                                  'models.pxi', 'fast_choose.pxi'],
                         language='c++',
                         include_dirs=[numpy.get_include()],
                         extra_compile_args=openmp_args,
                         extra_link_args=openmp_args,
    )

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

