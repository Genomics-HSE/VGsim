from setuptools import setup, Extension
import sys

boost_include = '/opt/homebrew/include'
boost_library = '/opt/homebrew/lib'
python_version = sys.version_info
boost_version = f'boost_python{python_version.major}{python_version.minor}'

example_module = Extension('source_VGsim',
                            sources=['VGsim_new/wrap.cpp'],
                            include_dirs=[boost_include],
                            library_dirs=[boost_library],
                            libraries=[boost_version],
                            extra_compile_args=['-std=c++17', '-O2'],
                            )

setup(name='VGsim_new',
      version='1.0',
      description='VGsim is the fast viral genealogy simulator for world-wide pandemic scenarios.',
      ext_modules=[example_module],
      packages=['VGsim_new'],
      install_requires=[
            'tabulate',
            'pytest',
            'numpy',
            'matplotlib'
      ])
