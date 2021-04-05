# ViralSim


Building the package
--------------------

For performance reasons, the package includes C extensions, which need to be
built. In the future, we plan to provide prebuilt packages, but for now you need
to build them yourself. For this, you need a working toolchain for building C++
code (gcc and clang are known to work). 

You can use either `pip` or `conda`. With `pip`, proceed as follows
(with `conda`, the process should be similar):
First, install the dependencies

```
$ python -m pip install numpy>=1.19.5 cython
$ python -m pip install git+https://github.com/ev-br/mc_lib.git@v0.1
```

Then, build the C extensions,

```
$ python setup.py build_ext --inplace
```

That's it! You may now run simulations:

```
$ python ./ViralSim.py example/example.rt -it 100000 -pm example/example.pp example/example.mg -seed 2020
```

If you encounter problems with either of these steps, please file an issue at
`https://github.com/Genomics-HSE/ViralSim` and include the build log.


We tested this procedure on python 3.7-3.9 on Ubuntu linux and MacOS. Whether
it works on Apple Silicon hardware, we do not know (most likely, it should
as soon as there is a NumPy version which supports this hardware).

![linux tests](https://github.com/ev-br/mc_lib/actions/workflows/python-package.yml/badge.svg)
![macOS tests](https://github.com/ev-br/mc_lib/actions/workflows/macos.yml/badge.svg)


The rest of the README
----------------------

File with rates is required.

Example:
```
python setup.py build_ext --inplace
./ViralSim.py example/example.rt -it 100000000 -pm example/example.pp example/example.mg -su example/example.su

for n in {1..5}; do ./ViralSim.py example/example.rt -it 1000000 -pm example/example.pp example/example.mg -su example/example.su; done >> errors.txt
```
