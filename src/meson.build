### Build the cython extension ###

# compile .cpp for _BirthDeath.pyx
_birth_death = custom_target('_BirthDeath',
  output : '_BirthDeath.cxx',
  input : '_BirthDeath.pyx',
  depend_files: ['models.pxi', 'fast_choose.pxi'],
  command : [cython, '--cplus', '-3', '--fast-fail', '@INPUT@', '-o', '@OUTPUT@']
)

npymath_path = incdir_numpy / '..' / 'lib'
npymath_lib = meson.get_compiler('c').find_library('npymath', dirs: npymath_path)
npyrandom_path = incdir_numpy / '..' / '..' / 'random' / 'lib'
npyrandom_lib = meson.get_compiler('c').find_library('npyrandom', dirs: npyrandom_path)

# generate .so for rndm
py3.extension_module(
  '_BirthDeath', _birth_death,
  include_directories: [inc_np, npymath_path, npyrandom_path],
  dependencies : [py3_dep, npymath_lib, npyrandom_lib],
  install: true,
  subdir: 'VGsim'
)


### Install stuff: python files + the extension just built ###

#create a list of files which will be installed as source
python_sources = [
    '__init__.py',
    'IO.py',
    '_interface.py',
]

# install source files in current dir
py3.install_sources(
    python_sources,
    pure: false,
    subdir: 'VGsim'
)
