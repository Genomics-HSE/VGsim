Various Examples
================

```{code-block} none

import VGsim.interface as inter

a = inter.Simulator()

a.initialize(2020)
a.simulate(200000)

```

```{code-block} none

$ conda create --name msprime-env
...

  environment location: /home/jean/miniconda3/envs/msprime-env

Proceed ([y]/n)? y

Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done
#
# To activate this environment, use
#
#     $ conda activate msprime-env
#
# To deactivate an active environment, use
#
#     $ conda deactivate

$ conda activate msprime-env
(msprime-env) $ conda install -c conda-forge msprime
...

  added / updated specs:
    - msprime

The following NEW packages will be INSTALLED:
...

Proceed ([y]/n)? y

Downloading and Extracting Packages

Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done

(msprime-env) $ python
Python 3.8.5
[GCC 7.3.0] :: Anaconda, Inc. on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import msprime
>>> ts = msprime.sim_ancestry(3)
>>> print(ts.draw_text())
3.18┊      10     ┊  
    ┊    ┏━━┻━━┓  ┊  
2.08┊    9     ┃  ┊  
    ┊  ┏━┻━┓   ┃  ┊  
0.80┊  ┃   ┃   8  ┊  
    ┊  ┃   ┃  ┏┻┓ ┊  
0.42┊  ┃   7  ┃ ┃ ┊  
    ┊  ┃  ┏┻┓ ┃ ┃ ┊  
0.34┊  6  ┃ ┃ ┃ ┃ ┊  
    ┊ ┏┻┓ ┃ ┃ ┃ ┃ ┊  
0.00┊ 0 4 1 3 2 5 ┊  
  0.00          1.00 
```