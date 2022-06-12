Building the package
--------------------

For performance reasons, the package includes C extensions, which need to be
built. In the future, we plan to provide prebuilt packages, but for now you need
to build them yourself. For this, you need a working toolchain for building C++
code (gcc and clang are known to work). Since you are going to build Python extensions,
you will need python development headers (e.g. on ubuntu linux the package name is `python-dev`).

The simplest quick-start cross-platform way is to use `conda`. To do this, create a fresh conda environment:

.. code-block:: bash

	$ conda create -n conda_vgsim
	$ conda activate conda_vgsim
	$ conda install python=3.9       # or other python version of your choice (any `python >= 3.7` should work).

Note that we install python *with conda, inside the conda environment*. Mixing system-install python and conda may lead to build- or runtime errors. 

Get the source code --- here we clone it from GitHub

.. code-block:: bash

	$ git clone https://github.com/Genomics-HSE/VGsim.git
	$ cd VGsim

Then build the package: 

.. code-block:: bash

	$ python -m pip install .

That's it! 
You may now run simulations:

.. code-block:: bash

	$ python ./testing/example.py

Running this script should create a directory ``testing/example_output`` with five files ``tree.nwk`` (genealogy tree), ``mutations.tsv`` (mutations on the genealogy), ``migrations.tsv`` (migrations of the lineages in the genealogy), ``sample_population.tsv`` (population of each sample), ``plot.png`` (visualization of epidemiological curves). There will be a warning in the output. This is an expected behaviour in this example (`see for detailed explanation <https://vg-sim.readthedocs.io/en/latest/Migration.html>`_ of inflated number of individuals in a deme due to migration). This example is also provided as a jupyter notebook ``testing/example.ipynb``.

If you prefer to not use ``conda``, the package builds fine with just pip. We recommend using virtual environments, via the standard library ``venv`` package (or a third-party ``virtualenv`` package). For MacOS users: please make sure to use Python with Homebrew. Avoid system Python. Manually installed Python (e.g. downloaded from python.org) currently does not work with meson too. If there are still errors, try to install pkg-config with brew install. YMMV.
We tested this procedure on python 3.7-3.9 on Ubuntu linux and MacOS.
On Apple Silicon, you need to have ``numpy >= 1.21`` (which is the first NumPy
version to support this hardware).

If you encounter problems with either of these steps, please file an issue at
``https://github.com/Genomics-HSE/VGsim``: please rerun with the ``-v`` flag,
``$ python -mpip install . -v`` and include the output.
