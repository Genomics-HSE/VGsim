from setuptools import setup, Extension
import sysconfig

# Найти путь к Boost
boost_include = '/usr/local/include'  # или путь, где установлен Boost
boost_library = '/usr/local/lib'  # или путь, где установлен Boost

example_module = Extension('VGsim_test', sources=['wrap.cpp'],
                            include_dirs=[boost_include],
                            library_dirs=[boost_library],
                            libraries=['boost_python312'],  # Убедитесь, что используете правильное имя библиотеки
                            extra_compile_args=['-std=c++17', '-O2'],
                            )

setup(name='VGsim_test',
      version='1.0',
      description='Example module built using Boost.Python',
      ext_modules=[example_module])
