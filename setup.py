from setuptools import setup, Extension
import sys

boost_include = '/usr/local/include'
boost_library = '/usr/local/lib'
python_version = sys.version_info
boost_version = f'boost_python{python_version.major}{python_version.minor}'

example_module = Extension('source_VGsim',
                            sources=['src/wrap.cpp'],
                            include_dirs=[boost_include],
                            library_dirs=[boost_library],
                            libraries=[boost_version],
                            extra_compile_args=['-std=c++17', '-O2', '-L/usr/local/lib', '-L/path/to/boost/libs', '-L/usr/local/include'],
                            )

setup(name='VGsim_new',
      version='1.0',
      description='VGsim is the fast viral genealogy simulator for world-wide pandemic scenarios.',
      ext_modules=[example_module],
      packages=['VGsim_new'])
