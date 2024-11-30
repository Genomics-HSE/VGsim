from setuptools import setup, Extension
import sysconfig
import sys

# Найти путь к Boost
# boost_include = '/usr/local/include'
# boost_library = '/usr/local/lib'
python_version = sys.version_info
boost_version = f'boost_python{python_version.major}{python_version.minor}'
# boost_version = f'boost_python{python_version.major}'


example_module = Extension('VGsim',
                            sources=['src/wrap.cpp'],
                            # include_dirs=[boost_include],
                            # library_dirs=[boost_library],
                            libraries=[boost_version],  # Убедитесь, что используете правильное имя библиотеки
                            extra_compile_args=['-std=c++17', '-O2'],
                            )

setup(name='VGsim',
      version='1.0',
      description='VGsim is the fast viral genealogy simulator for world-wide pandemic scenarios.',
      ext_modules=[example_module])
