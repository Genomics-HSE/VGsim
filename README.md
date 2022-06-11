# VGsim
Fast simulator of viral genealogies in the world-scale pandemic scenarios.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![pip](https://github.com/Genomics-HSE/VGsim/actions/workflows/pip.yml/badge.svg)](https://github.com/Genomics-HSE/VGsim/actions/workflows/pip.yml)
[![dev](https://github.com/Genomics-HSE/VGsim/actions/workflows/dev.yml/badge.svg)](https://github.com/Genomics-HSE/VGsim/actions/workflows/dev.yml)
[![windows](https://github.com/Genomics-HSE/VGsim/actions/workflows/windows.yml/badge.svg)](https://github.com/Genomics-HSE/VGsim/actions/workflows/windows.yml)
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

The simplest quick-start cross-platform way is to use `conda`. To do this, create a fresh conda environment:

```
$ conda create -n conda_vgsim
$ conda activate conda_vgsim
$ conda install python=3.9       # or other python version of your choice (any `python >= 3.7` should work).
```

Note that we install python *with conda, inside the conda environment*. Mixing system-install python and conda may lead to build- or runtime errors. 

Get the source code --- here we clone it from GitHub

```
$ git clone https://github.com/Genomics-HSE/VGsim.git
$ cd VGsim
```

Then build the package: 

```
$ python -m pip install .
```

That's it! 
You may now run simulations:

```
$ python ./testing/example.py
```

Running this script should create a directory ``testing/example_output`` with five files ``tree.nwk`` (genealogy tree), ``mutations.tsv`` (mutations on the genealogy), ``migrations.tsv`` (migrations of the lineages in the genealogy), ``sample_population.tsv`` (population of each sample), ``plot.png`` (visualization of epidemiological curves). There will be a warning in the output. This is an expected behaviour in this example (see https://vg-sim.readthedocs.io/en/latest/Migration.html for detailed explanation of inflated number of individuals in a deme due to migration). This example is also provided as a jupyter notebook ``testing/example.ipynb``.

If you prefer not to use ``conda``, the package builds fine with just pip. We recommend using virtual environments, via the standard library ``venv`` package (or a third-party ``virtualenv`` package). For MacOS users: please make sure to use Python with Homebrew. Avoid system Python. Manually installed Python (e.g. downloaded from python.org) currently does not work with meson too. If there are still errors, try to install pkg-config with brew install. YMMV.
We tested this procedure on python 3.7-3.9 on Ubuntu linux and MacOS.
On Apple Silicon, you need to have ``numpy >= 1.21`` (which is the first NumPy
version to support this hardware).

If you encounter problems with either of these steps, please file an issue at
``https://github.com/Genomics-HSE/VGsim``: please rerun with the ``-v`` flag,
``$ python -mpip install . -v`` and include the output.



Adding neutral mutations
------------------------

We suggest to pipe the tree and non-neautral sites into phastSim (https://github.com/NicolaDM/phastSim) by Nicola de Maio to add neutral mutations and to obtain full sequences.

Planned features
----------------
Here we list the features which we are adding or consider for the VGsim package.
- Tau-leaping: done (beta-version).
- Memory usage optimisation to enable larger number of mutable sites: done (beta-version).
- Population-level susceptibility transitions (e.g. to simulate different uneven efforts across the world): in work.
- Advanced immunity model: in work.
- Advanced sampling schemes: medium priority.
- Super-spreading events: medium priority.
- Life-cycle: low priority.
- Recombinations: low priority.


Citation
--------

Please cite our preprint (https://www.medrxiv.org/content/10.1101/2021.04.21.21255891) when using VGsim in your research.
```
Shchur et al. VGsim: scalable viral genealogy simulator for global pandemic.
medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891
```

Fast build
----------

For developers, it is recommended to manually build the package. This process is much faster then using pip, so it is handy while actively working with the code.

From the `VGsim` directory, first let us configure `PYTHOPATH`.
```
ver=$(python -c'from sys import version_info as v; print("%s.%s"%(v.major, v.minor))')
export PYTHONPATH=$PWD/installdir/lib/python$ver/site-packages/  
```
Remember to do it in every terminal window/tab where you work on the project (e.g. with your `jupyter notebook`)! Now we are ready to build `VGsim`.
```
meson setup build --prefix=$PWD/installdir    
meson install -C build
```
That's it! For subsequent rebuilds (in the same terminal session) you need to run only the last command
```
meson install -C build
```
