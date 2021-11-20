.. VGsim documentation master file, created by
   sphinx-quickstart on Wed Sep 29 16:30:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to VGsim's documentation!
=================================

VGsim is the Python package for simulation of viral genealogies in the world-scale pandemic scenarios designed primarily for research purposes. Key VGsim features include

* scalability and efficiency: optimized Guillespe algorithm over compartmental model generates epidemiological process; genealogy (transmission chain) is simulated only for sampled cases;
* pathogen evolution and epistasis (through single nucleotide substitutions): each haplotype properties (e.g. transmission rate, recovery rate, triggered host immunity) can be set independently;
* host immunity: different immunity types can have lower or higher susceptilibities to different haplotypes;
* population model: demic model with migration (migration is modelled as short-term travels);
* contact density (to reflect social, cultural and other population differences) and lockdowns.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   modules 

Citation
--------

Please cite our preprint (https://www.medrxiv.org/content/10.1101/2021.04.21.21255891) when using VGsim in your research.

Shchur et al. VGsim: scalable viral genealogy simulator for global pandemic.
Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig
medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891

