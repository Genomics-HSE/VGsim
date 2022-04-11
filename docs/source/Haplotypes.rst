Haplotypes
==========

One of the key features of VGsim is a usage of pathogen haplotypes (or strains). Haplotype determines haplotype properties: transmissibility, mutation rates, interaction with immunity etc. Currently VGsim is optimized for a relatively small number of haplotypes within a simulation. This number is equal to 4^U, where U is the number of genetic sites. Each site is assumed to have a single nucleotide with four possible alleles (A, T, C or G).

.. note::
    If you have many haplotypes in your simulation, we encourage you to set *memory_optimisation* field in **Simulator** class to *True*. That might help with the memory performance for your simulation.


Single nucleotide mutations lead to different haplotypes arising in the population. The rates of such mutations can be finely tuned: the user can set them independently for each site and each derived state of each haplotype.

VGsim’s API provides a convenient tool to set values of different parameters to a particular haplotype or to groups of haplotypes. Every API function with a parameter ``haplotype`` supports the following parameter values.

``haplotype`` default value is None: the setting will be applied to all possible haplotypes in the simulation.

First possibility to refer a haplotype or a group of haplotypes is to use a string of length U (U is the number of genetic sites in the simulation) over five possible symbols ‘ATCG*’. Here * asterisk corresponds to any nucleotide. For example, ‘*A*’ corresponds to 16 different haplotypes ‘AAA’, ‘AAT’, ‘AAC’, ‘AAG’, ‘TAA’, ‘TAT’, ‘TAC’, ‘TAG’ etc.

Another possibility is to use haplotype id, which is an integer. It can be obtained by converting a number in the quaternary numeral system (we assume the following bijection ATCG->0123) into the decimal system. In other words, haplotype CAG corresponds to a quaternary number 203. In decimal system it is 2*4^2+0^4^1+3*4^0=32+0+3=35.

Example
-------
Let us consider a handy way to set mutation rates. In this example we change the parameters for the mutation arising at the second site of four haplotypes AAA, ATA, ACA, AGA. The mutation rate is set to 2.5 for all cases. 

set_mutation_rate(2.5, [2, 2, 3, 1], haplotype=”A*A”, site_id=1)

The substitution weights will be the following:
	
Haplotype AAA: A **A** A->A **T** A: 2, A **A** A->A **C** A: 3, A **A** A->A **G** A:1

Haplotype ATA: A **T** A->A **A** A: 2, A **T** A->A **C** A: 3, A **T** A->A **G** A:1

Haplotype ACA: A **C** A->A **A** A: 2, A **C** A->A **T** A: 2, A **C** A->A **G** A:1

Haplotype AGA: A **G** A->A **A** A: 2, A **G** A->A **T** A: 3, A **G** A->A **C** A:3
