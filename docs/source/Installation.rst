Building the package
--------------------

For performance reasons, the package includes C extensions, which need to be
built. In the future, we plan to provide prebuilt packages, but for now you need
to build them yourself. For this, you need a working toolchain for building C++
code (gcc and clang are known to work). Since you are going to build Python extensions,
you will need python development headers (e.g. on ubuntu linux the package name is `python-dev`).


To build the C extensions, run

``$ python -m pip install .``

That's it! 

You may now run simulations:

``$ python ./vgsim.py example/example.rt -it 100000 -pm example/example.pp example/example.mg -seed 2020``

If you encounter problems with either of these steps, please file an issue at
``https://github.com/Genomics-HSE/VGsim``: please rerun with the ``-v`` flag,
``$ python -mpip install . -v`` and include the output.

We tested this procedure on python 3.7-3.9 on Ubuntu linux and MacOS. 
On Apple Silicon, you need to have `numpy >= 1.21` (which is the first NumPy
version to support this hardware).s