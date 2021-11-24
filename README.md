# VGsim
Fast simulator of viral genealogies in the world-scale pandemic scenarios.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pip](https://github.com/Genomics-HSE/VGsim/actions/workflows/pip.yml/badge.svg)](https://github.com/Genomics-HSE/VGsim/actions/workflows/pip.yml)
[![Documentation Status](https://readthedocs.org/projects/vg-sim/badge/?version=latest)](https://vg-sim.readthedocs.io/en/latest/?badge=latest)

Documentation
-------------

We keep on updating the detailed documentation (https://vg-sim.readthedocs.io/). In particular, there is a Python API which allows to build complicated simulations, which are currently not supported from the command line.

Building the package
--------------------

For performance reasons, the package includes C extensions, which need to be
built. In the future, we plan to provide prebuilt packages, but for now you need
to build them yourself. For this, you need a working toolchain for building C++
code (gcc and clang are known to work). Since you are going to build Python extensions,
you will need python development headers (e.g. on ubuntu linux the package name is `python-dev`).

**For MacOS users: please make sure to use Python with Homebrew or Conda. Avoid system Python. Manually installed Python (e.g. downloaded from python.org) currently does not work with meson too. If there are still errors, try install pkg-config with brew install**

To build the C extensions, run

```
$ python -m pip install .
```

That's it!

You may now run simulations:

```
$ python ./vgsim.py example/example.rt -it 100000 -pm example/example.pp example/example.mg -seed 2020
```

If you encounter problems with either of these steps, please file an issue at
`https://github.com/Genomics-HSE/VGsim`: please rerun with the `-v` flag,
`$ python -mpip install . -v` and include the output.

We tested this procedure on python 3.7-3.9 on Ubuntu linux and MacOS.
On Apple Silicon, you need to have `numpy >= 1.21` (which is the first NumPy
version to support this hardware).

Adding neutral mutations
------------------------

We suggest to pipe the tree and non-neautral sites into phastSim (https://github.com/NicolaDM/phastSim) by Nicola de Maio to add neutral mutations and to obtain full sequences.

Citation
--------

Please cite our preprint (https://www.medrxiv.org/content/10.1101/2021.04.21.21255891) when using VGsim in your research.
```
Shchur et al. VGsim: scalable viral genealogy simulator for global pandemic.
medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891
```
