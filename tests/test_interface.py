from numpy.testing import assert_allclose
import numpy as np

from VGsim import Simulator

import pytest


#SIMULATOR
@pytest.mark.parametrize('number_of_sites, answer', [(2, 2),
													 (0, 0)])
def test_simulator_sites(number_of_sites, answer):
	model = Simulator(number_of_sites=number_of_sites)
	assert model.number_of_sites == answer

@pytest.mark.parametrize('number_of_sites, error, text', [(None , TypeError , 'Incorrect type of number of sites. Type should be int.'),
														  ('str', TypeError , 'Incorrect type of number of sites. Type should be int.'),
														  (-2   , ValueError, 'Incorrect value of number of sites. Value should be more or equal 0.')])
def test_simulator_sites_err(number_of_sites, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(number_of_sites=number_of_sites)

@pytest.mark.parametrize('populations_number, answer', [(2, 2)])
def test_simulator_populations(populations_number, answer):
	model = Simulator(populations_number=populations_number)
	assert model.populations_number == answer

@pytest.mark.parametrize('populations_number, error, text', [(None , TypeError , 'Incorrect type of populations number. Type should be int.'),
															 ('str', TypeError , 'Incorrect type of populations number. Type should be int.'),
															 (0    , ValueError, 'Incorrect value of populations number. Value should be more 0.'),
															 (-2   , ValueError, 'Incorrect value of populations number. Value should be more 0.')])
def test_simulator_populations_err(populations_number, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(populations_number=populations_number)

@pytest.mark.parametrize('number_of_susceptible_groups, answer', [(2, 2)])
def test_simulator_susceptible(number_of_susceptible_groups, answer):
	model = Simulator(number_of_susceptible_groups=number_of_susceptible_groups)
	assert model.number_of_susceptible_groups == answer

@pytest.mark.parametrize('number_of_susceptible_groups, error, text', [(None , TypeError , 'Incorrect type of number of susceptible groups. Type should be int.'),
																	   ('str', TypeError , 'Incorrect type of number of susceptible groups. Type should be int.'),
																	   (0    , ValueError, 'Incorrect value of number of susceptible groups. Value should be more 0.'),
																	   (-2   , ValueError, 'Incorrect value of number of susceptible groups. Value should be more 0.')])
def test_simulator_susceptible_err(number_of_susceptible_groups, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(number_of_susceptible_groups=number_of_susceptible_groups)

@pytest.mark.parametrize('seed, answer', [(15, 15),
										  (0 , 0)])
def test_simulator_seed(seed, answer):
	model = Simulator(seed=seed)
	assert model.seed == answer

@pytest.mark.parametrize('seed, error, text', [('str', TypeError , 'Incorrect type of seed. Type should be int.'),
											  (-15  , ValueError, 'Incorrect value of seed. Value should be more or equal 0.')])
def test_simulator_seed_err(seed, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(seed=seed)


@pytest.mark.parametrize('sampling_probability, answer', [(True , True),
														  (False, False)])
def test_simulator_sampling_probability(sampling_probability, answer):
	model = Simulator(sampling_probability=sampling_probability)
	assert model.sampling_probability == answer

@pytest.mark.parametrize('sampling_probability, error, text', [(None , ValueError, 'Incorrect value of sampling probability. Value of sampling probability should be True or False.'),
															   ('str', ValueError, 'Incorrect value of sampling probability. Value of sampling probability should be True or False.')])
def test_simulator_sampling_probability_err(sampling_probability, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(sampling_probability=sampling_probability)

@pytest.mark.parametrize('memory_optimization, answer', [(True , True),
														 (False, False)])
def test_simulator_memory_optimization(memory_optimization, answer):
	model = Simulator(memory_optimization=memory_optimization)
	assert model.memory_optimization == answer

@pytest.mark.parametrize('memory_optimization, error, text', [(None , ValueError, 'Incorrect value of memory optimization. Value of memory optimization should be True or False.'),
															  ('str', ValueError, 'Incorrect value of memory optimization. Value of memory optimization should be True or False.')])
def test_simulator_memory_optimization_err(memory_optimization, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(memory_optimization=memory_optimization)

@pytest.mark.parametrize('genome_length, answer', [(1000, 1000)])
def test_simulator_genome_length(genome_length, answer):
	model = Simulator(genome_length=genome_length)
	assert model.genome_length == answer

@pytest.mark.parametrize('genome_length, error, text', [(None , TypeError , 'Incorrect type of genome length. Type should be int.'),
														('str', TypeError , 'Incorrect type of genome length. Type should be int.'),
														(-2   , ValueError, 'Incorrect value of genome length. Value should be more 0.'),
														(2    , ValueError, 'Incorrect value of number of sites or genome length. Genome length should be more or equal number of sites.')])
def test_simulator_genome_length_err(genome_length, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(number_of_sites=4, genome_length=genome_length)

@pytest.mark.parametrize('recombination_probability, answer', [(0  , 0),
															   (0.5, 0.5),
															   (1  , 1)])
def test_simulator_recombination_probability(recombination_probability, answer):
	model = Simulator(recombination_probability=recombination_probability)
	assert model.coinfection_parameters == answer


@pytest.mark.parametrize('recombination_probability, error, text', [(None , TypeError , 'Incorrect type of recombination probability. Type should be int or float.'),
																	('str', TypeError , 'Incorrect type of recombination probability. Type should be int or float.'),
																	(-1   , ValueError, 'Incorrect value of recombination probability. Value should be more or equal 0 and equal or less 1.'),
																	(1.1  , ValueError, 'Incorrect value of recombination probability. Value should be more or equal 0 and equal or less 1.')])
def test_simulator_recombination_probability_err(recombination_probability, error, text):
	with pytest.raises(error, match=text):
		model = Simulator(recombination_probability=recombination_probability)


#INITIAL HAPLOTYPE
@pytest.mark.parametrize('amount, answer', [(2 , 2),
											(18, 16)])
def test_set_initial_haplotype(amount, answer):
	model = Simulator(number_of_sites=2, memory_optimization=True)
	model.set_initial_haplotype(amount)
	assert model.initial_haplotype == answer

@pytest.mark.parametrize('memory_optimization, amount, error, text', [(False, 2    , ValueError, r'Incorrect value of memory optimization. Value should be equal \'True\' for work this function.'),
																	  (True , None , TypeError , 'Incorrect type of amount of initial haplotype. Type should be int.'),
																	  (True , 'str', TypeError , 'Incorrect type of amount of initial haplotype. Type should be int.'),
																	  (True , 0    , ValueError, 'Incorrect value of amount of initial haplotype. Value should be more 0.')])

def test_set_initial_haplotype_err(memory_optimization, amount, error, text):
	model = Simulator(number_of_sites=2, memory_optimization=memory_optimization)
	with pytest.raises(error, match=text):
		model.set_initial_haplotype(amount)


#STEP HAPLOTYPE
@pytest.mark.parametrize('amount, answer', [(2 , 2),
											(18, 18)])
def test_set_step_haplotype(amount, answer):
	model = Simulator(number_of_sites=2, memory_optimization=True)
	model.set_step_haplotype(amount)
	assert model.step_haplotype == answer

@pytest.mark.parametrize('memory_optimization, amount, error, text', [(False, 2    , ValueError, r'Incorrect value of memory optimization. Value should be equal \'True\' for work this function.'),
																	  (True , None , TypeError , 'Incorrect type of amount of step haplotype. Type should be int.'),
																	  (True , 'str', TypeError , 'Incorrect type of amount of step haplotype. Type should be int.'),
																	  (True , 0    , ValueError, 'Incorrect value of amount of step haplotype. Value should be more 0.')])

def test_set_step_haplotype_err(memory_optimization, amount, error, text):
	model = Simulator(number_of_sites=2, memory_optimization=memory_optimization)
	with pytest.raises(error, match=text):
		model.set_step_haplotype(amount)


#GENOME LENGTH
@pytest.mark.parametrize('genome_length, answer', [(1000, 1000)])
def test_set_genome_length(genome_length, answer):
	model = Simulator()
	model.set_genome_length(genome_length=genome_length)
	assert model.genome_length == answer

@pytest.mark.parametrize('genome_length, error, text', [(None , TypeError , 'Incorrect type of genome length. Type should be int.'),
														('str', TypeError , 'Incorrect type of genome length. Type should be int.'),
														(-4   , ValueError, 'Incorrect value of genome length. Value should be more 0.'),
														(3    , ValueError, 'Incorrect value of number of sites or genome length. Genome length should be more or equal number of sites.')])
def test_set_genome_length_err(genome_length, error, text):
	model = Simulator(number_of_sites=5)
	with pytest.raises(error, match=text):
		model.set_genome_length(genome_length=genome_length)


#COINFECTION PARAMETERS
@pytest.mark.parametrize('recombination, answer', [(0.4, 0.4)])
def test_set_coinfection_parameters(recombination, answer):
	model = Simulator()
	model.set_coinfection_parameters(recombination=recombination)
	assert model.coinfection_parameters == answer

@pytest.mark.parametrize('recombination, error, text', [(None , TypeError , 'Incorrect type of recombination probability. Type should be int.'),
														('str', TypeError , 'Incorrect type of recombination probability. Type should be int.'),
														(-1   , ValueError, 'Incorrect value of recombination probability. Value should be more or equal 0 and equal or less 1.'),
														(1.2  , ValueError, 'Incorrect value of recombination probability. Value should be more or equal 0 and equal or less 1.')])
def test_set_coinfection_parameters(recombination, error, text):
	model = Simulator()
	with pytest.raises(error, match=text):
		model.set_coinfection_parameters(recombination=recombination)


#TRANSMISSION RATE
# @pytest.mark.parametrize('rate, haplotype, condition', [()])
@pytest.mark.parametrize('rate, haplotype, answer', [(0    , None, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
													 (0.002, 'A*', [0.002, 0.002, 0.002, 0.002, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]),
													 (0.003, 'AT', [2, 0.003, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]),
													 (0.004, 0   , [0.004, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]),
													 (0.005, 15  , [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0.005])])
def test_set_transmisstion_rate(rate, haplotype, answer):
	model = Simulator(number_of_sites=2)
	model.set_transmission_rate(rate=rate, haplotype=haplotype)
	# print(type(model.get_indexes_from_haplotype(haplotype)[0]))
	# model.transmission_rate[model.get_indexes_from_haplotype(haplotype)] = rate  #(rate=rate, haplotype=haplotype)
	# assert_allclose(model.get_transmission_rate(), np.asarray(answer), atol=1e-14)
	assert_allclose(model.transmission_rate, np.asarray(answer), atol=1e-14)

# @pytest.mark.parametrize('rate, haplotype, condition, error', [()])
@pytest.mark.parametrize('rate, haplotype, error, text', [(None , None  , TypeError , 'Incorrect type of transmission rate. Type should be int or float.'),
														  ('str', None  , TypeError , 'Incorrect type of transmission rate. Type should be int or float.'),
														  (-1   , None  , ValueError, 'Incorrect value of transmission rate. Value should be more or equal 0.'),
														  (0.01 , [1, 2], TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
														  (0.01 , -1    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 16    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 'str' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  (0.01 , 'AAA' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.')]) 
def test_set_transmisstion_rate_err(rate, haplotype, error, text):
	model = Simulator(number_of_sites=2)
	with pytest.raises(error, match=text):
		model.set_transmission_rate(rate=rate, haplotype=haplotype)


#RECOVERY RATE
@pytest.mark.parametrize('rate, haplotype, answer', [(0    , None, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
													 (0.002, 'A*', [0.002, 0.002, 0.002, 0.002, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
													 (0.003, 'AT', [1, 0.003, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
													 (0.004, 0   , [0.004, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
													 (0.005, 15  , [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.005])])
def test_set_recovery_rate(rate, haplotype, answer):
	model = Simulator(number_of_sites=2)
	model.set_recovery_rate(rate=rate, haplotype=haplotype)
	assert_allclose(model.recovery_rate, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('rate, haplotype, error, text', [(None , None  , TypeError , 'Incorrect type of recovery rate. Type should be int or float.'),
														  ('str', None  , TypeError , 'Incorrect type of recovery rate. Type should be int or float.'),
														  (-1   , None  , ValueError, 'Incorrect value of recovery rate. Value should be more or equal 0.'),
														  (0.01 , [1, 2], TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
														  (0.01 , -1    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 16    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 'str' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  (0.01 , 'AAA' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.')]) 
def test_set_recovery_rate_err(rate, haplotype, error, text):
	model = Simulator(number_of_sites=2)
	with pytest.raises(error, match=text):
		model.set_recovery_rate(rate=rate, haplotype=haplotype)

#SAMPLING RATE
@pytest.mark.parametrize('rate, haplotype, answer', [(0    , None, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
													 (0.002, 'A*', [0.002, 0.002, 0.002, 0.002, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]),
													 (0.003, 'AT', [0.01, 0.003, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]),
													 (0.004, 0   , [0.004, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]),
													 (0.005, 15  , [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.005])])
def test_set_sampling_rate(rate, haplotype, answer):
	model = Simulator(number_of_sites=2)
	model.set_sampling_rate(rate=rate, haplotype=haplotype)
	assert_allclose(model.sampling_rate, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('rate, haplotype, error, text', [(None , None  , TypeError , 'Incorrect type of sampling rate. Type should be int or float.'),
														  ('str', None  , TypeError , 'Incorrect type of sampling rate. Type should be int or float.'),
														  (-1   , None  , ValueError, 'Incorrect value of sampling rate. Value should be more or equal 0.'),
														  (0.01 , [1, 2], TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
														  (0.01 , -1    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 16    , IndexError, 'There are no such haplotype!'),
														  (0.01 , 'str' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  (0.01 , 'AAA' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.')]) 
def test_set_sampling_rate_err(rate, haplotype, error, text):
	model = Simulator(number_of_sites=2)
	with pytest.raises(error, match=text):
		model.set_sampling_rate(rate=rate, haplotype=haplotype)

# #MUTATION RATE
@pytest.mark.parametrize('rate, haplotype, mutation, answer', [(0.001, None, None, [[0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001], [0.001, 0.001]]),
															   (0.002, 'A*', None, [[0.002, 0.002], [0.002, 0.002], [0.002, 0.002], [0.002, 0.002], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.003, 'AT', None, [[0.01, 0.01], [0.003, 0.003], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.004, 0   , None, [[0.004, 0.004], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.005, 15  , None, [[0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.005, 0.005]]),
															   (0.006, None, 0   , [[0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01], [0.006, 0.01]]),
															   (0.007, 'A*', 0   , [[0.007, 0.01], [0.007, 0.01], [0.007, 0.01], [0.007, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.008, 'AT', 0   , [[0.01, 0.01], [0.008, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.009, 0   , 0   , [[0.009, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01]]),
															   (0.011, 15  , 0   , [[0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.01, 0.01], [0.011, 0.01]])])
def test_set_mutation_rate(rate, haplotype, mutation, answer):
	model = Simulator(number_of_sites=2)
	model.set_mutation_rate(rate=rate, haplotype=haplotype, mutation=mutation)
	assert_allclose(model.mutation_rate, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('rate, haplotype, mutation, error, text', [(None , None  , None , TypeError , 'Incorrect type of mutation rate. Type should be int or float.'),
																	('str', None  , None , TypeError , 'Incorrect type of mutation rate. Type should be int or float.'),
																	(-1   , None  , None , ValueError, 'Incorrect value of mutation rate. Value should be more or equal 0.'),
																	(2    , [1, 2], None , TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
																	(2    , -1    , None , IndexError, 'There are no such haplotype!'),
																	(2    , 16    , None , IndexError, 'There are no such haplotype!'),
																	(2    , 'str' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  			(2    , 'AAA' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  			(2    , None  , 'str', TypeError , 'Incorrect type of mutation site. Type should be int or None.'),
														  			(2    , None  , -1   , IndexError, 'There are no such mutation site!'),
														  			(2    , None  , 2    , IndexError, 'There are no such mutation site!')])
def test_set_mutation_rate(rate, haplotype, mutation, error, text):
	model = Simulator(number_of_sites=2)
	with pytest.raises(error, match=text):
		model.set_mutation_rate(rate=rate, haplotype=haplotype, mutation=mutation)


#MUTATION PROBABILITY
@pytest.mark.parametrize('probabilities, haplotype, mutation, answer', [([2, 3, 4, 5], None, None, [[[3.0, 4.0, 5.0], [3.0, 4.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 4.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 3.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 3.0, 4.0]],[[2.0, 4.0, 5.0], [3.0, 4.0, 5.0]],[[2.0, 4.0, 5.0], [2.0, 4.0, 5.0]],[[2.0, 4.0, 5.0], [2.0, 3.0, 5.0]],[[2.0, 4.0, 5.0], [2.0, 3.0, 4.0]],[[2.0, 3.0, 5.0], [3.0, 4.0, 5.0]],[[2.0, 3.0, 5.0], [2.0, 4.0, 5.0]],[[2.0, 3.0, 5.0], [2.0, 3.0, 5.0]],[[2.0, 3.0, 5.0], [2.0, 3.0, 4.0]],[[2.0, 3.0, 4.0], [3.0, 4.0, 5.0]],[[2.0, 3.0, 4.0], [2.0, 4.0, 5.0]],[[2.0, 3.0, 4.0], [2.0, 3.0, 5.0]],[[2.0, 3.0, 4.0], [2.0, 3.0, 4.0]]]),
																		([2, 3, 4, 5], 'A*', None, [[[3.0, 4.0, 5.0], [3.0, 4.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 4.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 3.0, 5.0]],[[3.0, 4.0, 5.0], [2.0, 3.0, 4.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 'AT', None, [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [2.0, 4.0, 5.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 0   , None, [[[3.0, 4.0, 5.0], [3.0, 4.0, 5.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 15  , None, [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [2.0, 3.0, 4.0]]]),
																		([2, 3, 4, 5], None, 0   , [[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 5.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 'A*', 0   , [[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 'AT', 0   , [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 0   , 0   , [[[3.0, 4.0, 5.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]]),
																		([2, 3, 4, 5], 15  , 0   , [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],[[2.0, 3.0, 4.0], [1.0, 1.0, 1.0]]])])
def test_set_mutation_probabilities(probabilities, haplotype, mutation, answer):
	model = Simulator(number_of_sites=2)
	model.set_mutation_probabilities(probabilities=probabilities, haplotype=haplotype, mutation=mutation)
	assert_allclose(model.mutation_probabilities, np.asarray(answer), atol=1e-14)


@pytest.mark.parametrize('probabilities, haplotype, mutation, error, text', [(None            , None  , None , TypeError , 'Incorrect type of probabilities list. Type should be list.'),
																			 ('str'           , None  , None , TypeError , 'Incorrect type of probabilities list. Type should be list.'),
																			 ([0, 0, 0]       , None  , None , ValueError, 'Incorrect length of probabilities list. Length should be equal 4.'),
																			 ([0, 0, 0, 5]    , None  , None , ValueError, 'Incorrect probabilities list. The sum of three elements without mutation allele should be more 0.'),
																			 ([None, 0, 0, 5] , None  , None , TypeError , 'Incorrect type of mutation probabilities. Type should be int or float.'),
																			 (['str', 0, 0, 5], None  , None , TypeError , 'Incorrect type of mutation probabilities. Type should be int or float.'),
																			 ([2, 3, 4, 5]    , [1, 2], None , TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
																			 ([2, 3, 4, 5]    , -1    , None , IndexError, 'There are no such haplotype!'),
																			 ([2, 3, 4, 5]    , 16    , None , IndexError, 'There are no such haplotype!'),
																			 ([2, 3, 4, 5]    , 'str' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  					 ([2, 3, 4, 5]    , 'AAA' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
														  					 ([2, 3, 4, 5]    , None  , 'str', TypeError , 'Incorrect type of mutation site. Type should be int or None.'),
														  					 ([2, 3, 4, 5]    , None  , -1   , IndexError, 'There are no such mutation site!'),
														  					 ([2, 3, 4, 5]    , None  , 2    , IndexError, 'There are no such mutation site!')])
def test_set_mutation_probabilities_err(probabilities, haplotype, mutation, error, text):
	model = Simulator(number_of_sites=2)
	with pytest.raises(error, match=text):
		model.set_mutation_probabilities(probabilities=probabilities, haplotype=haplotype, mutation=mutation)


#MUTATION POSITION
@pytest.mark.parametrize('mutation, position, answer', [(2, 1000, [0, 333333, 1000, 1000000])])
def test_set_mutation_position(mutation, position, answer):
	model = Simulator(number_of_sites=4)
	model.set_mutation_position(mutation=mutation, position=position)
	assert_allclose(model.mutation_position, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('mutation, position, error, text', [(None , 1000   , TypeError , 'Incorrect type of number of site. Type should be int.'),
															 ('str', 1000   , TypeError , 'Incorrect type of number of site. Type should be int.'),
															 (-1   , 1000   , IndexError, 'There are no such number of site!'),
															 (4    , 1000   , IndexError, 'There are no such number of site!'),
															 (2    , None   , TypeError , 'Incorrect type of mutation position. Type should be int.'),
															 (2    , 'str'  , TypeError , 'Incorrect type of mutation position. Type should be int.'),
															 (2    , -1     , IndexError, 'There are no such mutation position!'),
															 (2    , 1000001, IndexError, 'There are no such mutation position!'),
															 (2    , 0      , IndexError, 'Incorrect value of position. Two mutations can\'t have the same position.')])
def test_set_mutation_position_err(mutation, position, error, text):
	model = Simulator(number_of_sites=4)
	with pytest.raises(error, match=text):
		model.set_mutation_position(mutation=mutation, position=position)


#POPULATION SIZES
@pytest.mark.parametrize('size, population, answer', [(1000, None, [1000, 1000, 1000, 1000]),
													  (1001, 2   , [1000000, 1000000, 1001, 1000000])])
def test_set_population_size(size, population, answer):
	model = Simulator(populations_number=4)
	model.set_population_size(size=size, population=population)
	assert_allclose(model.population_size, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('size, population, error, text', [(None , None , TypeError , 'Incorrect type of population size. Type should be int.'),
														   (2.5  , None , TypeError , 'Incorrect type of population size. Type should be int.'),
														   ('str', None , TypeError , 'Incorrect type of population size. Type should be int.'),
														   (-1   , None , ValueError, 'Incorrect value of population size. Value should be more 0.'),
														   (1000 , 1.5  , TypeError , 'Incorrect type of population. Type should be int or None.'),
														   (1000 , 'str', TypeError , 'Incorrect type of population. Type should be int or None.'),
														   (1000 , -2   , IndexError, 'There are no such population!'),
														   (1000 , 5    , IndexError, 'There are no such population!')])
def test_set_population_size_err(size, population, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_population_size(size=size, population=population)


#CONTACT DENSITY
@pytest.mark.parametrize('value, population, answer', [(2  , None, [2, 2, 2, 2]),
													   (3  , 2   , [1, 1, 3, 1]),
													   (2.5, 3   , [1, 1, 1, 2.5])])
def test_set_contact_density(value, population, answer):
	model = Simulator(populations_number=4)
	model.set_contact_density(value=value, population=population)
	assert_allclose(model.contact_density, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('value, population, error, text', [(None , None , TypeError , 'Incorrect type of contact density. Type should be int.'),
														    ('str', None , TypeError , 'Incorrect type of contact density. Type should be int.'),
														    (-1   , None , ValueError, 'Incorrect value of contact density. Value should be more or equal 0.'),
														    (2    , 1.5  , TypeError , 'Incorrect type of population. Type should be int or None.'),
														    (2    , 'str', TypeError , 'Incorrect type of population. Type should be int or None.'),
														    (2    , -2   , IndexError, 'There are no such population!'),
														    (2    , 5    , IndexError, 'There are no such population!')])
def test_set_contact_density_err(value, population, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_contact_density(value=value, population=population)


#SAMPLING MULTIPLIER
@pytest.mark.parametrize('multiplier, population, answer', [(2  , None, [2, 2, 2, 2]),
															(3  , 2   , [1, 1, 3, 1]),
															(2.5, 3   , [1, 1, 1, 2.5])])
def test_set_sampling_multiplier(multiplier, population, answer):
	model = Simulator(populations_number=4)
	model.set_sampling_multiplier(multiplier=multiplier, population=population)
	assert_allclose(model.sampling_multiplier, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('multiplier, population, error, text', [(None , None , TypeError , 'Incorrect type of sampling multiplier. Type should be int.'),
																 ('str', None , TypeError , 'Incorrect type of sampling multiplier. Type should be int.'),
																 (-1   , None , ValueError, 'Incorrect value of sampling multiplier. Value should be more or equal 0.'),
																 (2    , 1.5  , TypeError , 'Incorrect type of population. Type should be int or None.'),
																 (2    , 'str', TypeError , 'Incorrect type of population. Type should be int or None.'),
																 (2    , -2   , IndexError, 'There are no such population!'),
																 (2    , 5    , IndexError, 'There are no such population!')])
def test_set_sampling_multiplier_err(multiplier, population, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_sampling_multiplier(multiplier=multiplier, population=population)

#NPI
@pytest.mark.parametrize('parameters, population, answer', [([2, 0.7, 0.3], None, [[2, 2, 2, 2], [0.7, 0.7, 0.7, 0.7], [0.3, 0.3, 0.3, 0.3]]),
															([2, 0.7, 0.3], 2   , [[0, 0, 2, 0], [1.0, 1.0, 0.7, 1.0], [1.0, 1.0, 0.3, 1.0]])])
def test_set_npi(parameters, population, answer):
	model = Simulator(populations_number=4)
	model.set_npi(parameters=parameters, population=population)
	assert_allclose(model.npi, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('parameters, population, error, text', [(None              , None , TypeError , 'Incorrect type of npi parameters. Type should be list.'),
																 ('str'             , None , TypeError , 'Incorrect type of npi parameters. Type should be list.'),
																 ([2, 0.5, 0.5, 0.5], None , ValueError, 'Incorrect length of npi parameters. Length should be equal 3.'),
																 ([None, 0.5, 0.5]  , None , TypeError , 'Incorrect type of first npi parameter. Type should be int or float.'),
																 (['str', 0.5, 0.5] , None , TypeError , 'Incorrect type of first npi parameter. Type should be int or float.'),
																 ([-2, 0.5, 0.5]    , None , ValueError, 'Incorrect value of first npi parameter. Value should be more or equal 0.'),
																 ([2, None, 0.5]    , None , TypeError , 'Incorrect type of second npi parameter. Type should be int or float.'),
																 ([2, 'str', 0.5]   , None , TypeError , 'Incorrect type of second npi parameter. Type should be int or float.'),
																 ([2, -0.2, 0.5]    , None , ValueError, 'Incorrect value of second npi parameter. Value should be more or equal 0 and equal or less 1.'),
																 ([2, 1.1, 0.5]     , None , ValueError, 'Incorrect value of second npi parameter. Value should be more or equal 0 and equal or less 1.'),
																 ([2, 0.5, None]    , None , TypeError , 'Incorrect type of third npi parameter. Type should be int or float.'),
																 ([2, 0.5, 'str']   , None , TypeError , 'Incorrect type of third npi parameter. Type should be int or float.'),
																 ([2, 0.5, -0.2]    , None , ValueError, 'Incorrect value of third npi parameter. Value should be more or equal 0 and equal or less 1.'),
																 ([2, 0.5, 1.2]     , None , ValueError, 'Incorrect value of third npi parameter. Value should be more or equal 0 and equal or less 1.'),
																 ([2, 0.3, 0.3]     , 'str', TypeError , 'Incorrect type of population. Type should be int or None.'),
																 ([2, 0.3, 0.3]     , -1   , IndexError, 'There are no such population!'),
																 ([2, 0.3, 0.3]     , 4    , IndexError, 'There are no such population!')])
def test_set_npi_err(parameters, population, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_npi(parameters=parameters, population=population)


#MIGRATION PROBABILITY
@pytest.mark.parametrize('probability, source, target, answer', [(0.1 , None, None, [[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.1, 0.1, 0.1, 0.7]]),
																 (0.2 , 2   , None, [[1, 0, 0, 0], [0, 1, 0, 0], [0.2, 0.2, 0.4, 0.2], [0, 0, 0, 1]]),
																 (0.3 , None, 3   , [[0.7, 0, 0, 0.3], [0, 0.7, 0, 0.3], [0, 0, 0.7, 0.3], [0, 0, 0, 1]]),
																 (0.25, 1   , 2   , [[1, 0, 0, 0], [0, 0.75, 0.25, 0], [0, 0, 1, 0], [0, 0, 0, 1]])])
def test_set_migration_probability(probability, source, target, answer):
	model = Simulator(populations_number=4)
	model.set_migration_probability(probability=probability, source=source, target=target)
	assert_allclose(model.migration_probability, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('probability, source, target, error, text', [(None   , None , None , TypeError  , 'Incorrect type of migration probability. Type should be int or float.'),
																	  ('str'  , None , None , TypeError  , 'Incorrect type of migration probability. Type should be int or float.'),
																	  (-0.1   , None , None , ValueError , 'Incorrect value of migration probability. Value should be more or equal 0 and equal or less 1.'),
																	  (0.4    , None , None , ValueError , 'Incorrect the sum of migration probabilities. The sum of migration probabilities from each population should be equal or less 1.'),
																	  (1.0/3  , None , None , ValueError , 'Incorrect value of migration probability. Value of migration probability from source population to target population should be more 0.'),
																	  (0.1    , 'str', None , TypeError  , 'Incorrect type of population. Type should be int or None.'),
																	  (0.1    , -1   , None , IndexError , 'There are no such population!'),
																	  (0.1    , 4    , None , IndexError , 'There are no such population!'),
																	  (0.1    , None , 'str', TypeError  , 'Incorrect type of population. Type should be int or None.'),
																	  (0.1    , None , -1   , IndexError , 'There are no such population!'),
																	  (0.1    , None , 4    , IndexError , 'There are no such population!')])
def test_set_migration_probability_err(probability, source, target, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_migration_probability(probability=probability, source=source, target=target)


#TOTAL MIGRATION PROBABILITY
@pytest.mark.parametrize('total_probability, answer', [(0.3, [[0.7, 0.1, 0.1, 0.1], [0.1, 0.7, 0.1, 0.1], [0.1, 0.1, 0.7, 0.1], [0.1, 0.1, 0.1, 0.7]])])
def test_set_total_migration_probability(total_probability, answer):
	model = Simulator(populations_number=4)
	model.set_total_migration_probability(total_probability=total_probability)
	assert_allclose(model.migration_probability, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('total_probability, error, text', [(None , TypeError , 'Incorrect type of total migration probability. Type should be int or float.'),
															('str', TypeError , 'Incorrect type of total migration probability. Type should be int or float.'),
															(-0.1 , ValueError, 'Incorrect value of total migration probability. Value should be more or equal 0 and equal or less 1.'),
															(1.2  , ValueError, 'Incorrect value of total migration probability. Value should be more or equal 0 and equal or less 1.'),
															(1.0  , ValueError, 'Incorrect value of migration probability. Value of migration probability from source population to target population should be more 0.')])
def test_set_total_migration_probability_err(total_probability, error, text):
	model = Simulator(populations_number=4)
	with pytest.raises(error, match=text):
		model.set_total_migration_probability(total_probability=total_probability)



#SUSCEPTIBILITY TYPE
@pytest.mark.parametrize('susceptibility_type, haplotype, answer', [(1, None, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
																	(1, 'A*', [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
																	(1, 'AT', [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
																	(1, 0   , [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
																	(1, 15  , [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])])
def test_set_susceptibility_type(susceptibility_type, haplotype, answer):
	model = Simulator(number_of_sites=2, number_of_susceptible_groups=4)
	model.set_susceptibility_type(susceptibility_type=susceptibility_type, haplotype=haplotype)
	assert_allclose(model.susceptibility_type, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('susceptibility_type, haplotype, error, text', [(None , None  , TypeError , 'Incorrect type of susceptibility type. Type should be int.'),
																		 ('str', None  , TypeError , 'Incorrect type of susceptibility type. Type should be int.'),
																		 (-1   , None  , IndexError, 'There are no such susceptibility type!'),
																		 (4    , None  , IndexError, 'There are no such susceptibility type!'),
																		 (2    , [1, 2], TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
																		 (2    , -1    , IndexError, 'There are no such haplotype!'),
																		 (2    , 16    , IndexError, 'There are no such haplotype!'),
																		 (2    , 'str' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
																		 (2    , 'AAA' , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.')])
def test_set_susceptibility_type_err(susceptibility_type, haplotype, error, text):
	model = Simulator(number_of_sites=2, number_of_susceptible_groups=4)
	with pytest.raises(error, match=text):
		model.set_susceptibility_type(susceptibility_type=susceptibility_type, haplotype=haplotype)


#SUSCEPTIBILITY
@pytest.mark.parametrize('rate, haplotype, susceptibility_type, answer', [(2, None, None, [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]]),
																		  (3, 'A*', None, [[3, 3, 3], [3, 3, 3], [3, 3, 3], [3, 3, 3], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]]),
																		  (4, 'AT', None, [[1, 0, 0], [4, 4, 4], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]]),
																		  (5, 0   , None, [[5, 5, 5], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]]),
																		  (6, 15  , None, [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [6, 6, 6]]),
																		  (7, None, 2   , [[1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7], [1, 0, 7]])])
def test_set_susceptibility(rate, haplotype, susceptibility_type, answer):
	model = Simulator(number_of_sites=2, number_of_susceptible_groups=3)
	model.set_susceptibility(rate=rate, haplotype=haplotype, susceptibility_type=susceptibility_type)
	assert_allclose(model.susceptibility, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('rate, haplotype, susceptibility_type, error, text', [(None , None  , None , TypeError , 'Incorrect type of susceptibility rate. Type should be int or float.'),
																			   ('str', None  , None , TypeError , 'Incorrect type of susceptibility rate. Type should be int or float.'),
																			   (-1   , None  , None , ValueError, 'Incorrect value of susceptibility rate. Value should be more or equal 0.'),
																			   (2    , [1, 2], None , TypeError , 'Incorrect type of haplotype. Type should be int or str or None.'),
																			   (2    , -1    , None , IndexError, 'There are no such haplotype!'),
																			   (2    , 16    , None , IndexError, 'There are no such haplotype!'),
																			   (2    , 'str' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
																			   (2    , 'AAA' , None , ValueError, r'Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"\*\" and length of haplotype should be equal number of mutations sites.'),
																			   (2    , None  , 'str', TypeError , 'Incorrect type of susceptibility type. Type should be int or None.'), 
																			   (2    , None  , -1   , IndexError, 'There are no such susceptibility type!'),
																			   (2    , None  , 4    , IndexError, 'There are no such susceptibility type!')])
def test_set_susceptibility_err(rate, haplotype, susceptibility_type, error, text):
	model = Simulator(number_of_sites=2, number_of_susceptible_groups=4)
	with pytest.raises(error, match=text):
		model.set_susceptibility(rate=rate, haplotype=haplotype, susceptibility_type=susceptibility_type)


#IMMUNITY TRANSITION
@pytest.mark.parametrize('rate, source, target, answer', [(1, None, None, [[0, 1, 1], [1, 0, 1], [1, 1, 0]]),
														  (2, 1   , None, [[0, 0, 0], [2, 0, 2], [0, 0, 0]]),
														  (3, None, 2   , [[0, 0, 3], [0, 0, 3], [0, 0, 0]]),
														  (4, 0   , 2   , [[0, 0, 4], [0, 0, 0], [0, 0, 0]])])
def test_set_immunity_transition(rate, source, target, answer):
	model = Simulator(number_of_susceptible_groups=3)
	model.set_immunity_transition(rate=rate, source=source, target=target)
	assert_allclose(model.immunity_transition, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('rate, source, target, error, text', [(None , None , None , TypeError , 'Incorrect type of immunity transition rate. Type should be int or float.'),
															   ('str', None , None , TypeError , 'Incorrect type of immunity transition rate. Type should be int or float.'),
															   (-1   , None , None , ValueError, 'Incorrect value of immunity transition rate. Value should be more or equal 0.'),
															   (2    , -2   , None , IndexError, 'There are no such susceptibility type!'),
															   (2    , 5    , None , IndexError, 'There are no such susceptibility type!'),
															   (2    , 1.5  , None , TypeError , 'Incorrect type of susceptibility type. Type should be int or None.'),
															   (2    , 'str', None , TypeError , 'Incorrect type of susceptibility type. Type should be int or None.'),
															   (2    , None , -2   , IndexError, 'There are no such susceptibility type!'),
															   (2    , None , 5    , IndexError, 'There are no such susceptibility type!'),
															   (2    , None , 1.5  , TypeError , 'Incorrect type of susceptibility type. Type should be int or None.'),
															   (2    , None , 'str', TypeError , 'Incorrect type of susceptibility type. Type should be int or None.')])
def test_set_immunity_transition_err(rate, source, target, error, text):
	model = Simulator(number_of_susceptible_groups=3)
	with pytest.raises(error, match=text):
		model.set_immunity_transition(rate=rate, source=source, target=target)


#GENERAL SAMPLING
@pytest.mark.parametrize('sampling_proportion, sampling_times, answer', [(0.1, [300, 200], [[200, 0.1], [300, 0.1]]),
													 (0.2, 100.5, [[100.5, 0.2]]),
													 (0, [400, 300, 200],[[200, 0], [300, 0], [400, 0]]),
													 (1, 500, [[500, 1]])])
def test_set_general_sampling(sampling_proportion, sampling_times, answer):
	model = Simulator()
	model.set_general_sampling(sampling_proportion=sampling_proportion, sampling_times=sampling_times)
	assert_allclose(model.general_sampling_conditions, np.asarray(answer), atol=1e-14)

@pytest.mark.parametrize('sampling_proportion, sampling_times, error, text', [(None, None, TypeError, 'Incorrect type of general sampling proportion. Type should be int or float.'),
														  ('str', None, TypeError , 'Incorrect type of general sampling proportion. Type should be int or float.'),
														  (-1, None, ValueError, 'Incorrect value of general sampling proportion. Value should be more or equal 0 and equal or less 1.'),
														  (2, None, ValueError, 'Incorrect value of general sampling proportion. Value should be more or equal 0 and equal or less 1.'),
														  (0.01, None, TypeError , 'Incorrect type of list of general sampling times. Type should be list, int or float.'),
														  (0.01, 'str', TypeError, 'Incorrect type of list of general sampling times. Type should be list, int or float.'),
														  (0.01, [], ValueError, 'Incorrect length of list of general sampling times. Length should be greater than 0.'),
														  (0.01, [None], TypeError, 'Incorrect type of general sampling time. Type should be int or float.'),
														  (0.01, ['str'], TypeError, 'Incorrect type of general sampling time. Type should be int or float.'),
														  (0.01, -100, ValueError, 'Incorrect value of general sampling time. Value should be more or equal 0.')])
def test_set_general_sampling_err(sampling_proportion, sampling_times, error, text):
	model = Simulator()
	with pytest.raises(error, match=text):
		model.set_general_sampling(sampling_proportion=sampling_proportion, sampling_times=sampling_times)