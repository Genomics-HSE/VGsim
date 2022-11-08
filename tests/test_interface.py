from numpy.testing import assert_allclose
import numpy as np

from VGsim import Simulator

import pytest

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
	# assert_allclose(model.get_transmission_rate(), np.asarray(answer), atol=1e-14)
	assert_allclose(model.transmission_rate, np.asarray(answer), atol=1e-14)

# @pytest.mark.parametrize('rate, haplotype, condition, error', [()])
@pytest.mark.parametrize('rate, haplotype, error, text', [('str', None  , TypeError , 'Incorrect type of transmission rate. Type should be int or float.'),
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

@pytest.mark.parametrize('rate, haplotype, error, text', [('str', None  , TypeError , 'Incorrect type of recovery rate. Type should be int or float.'),
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

@pytest.mark.parametrize('rate, haplotype, error, text', [('str', None  , TypeError , 'Incorrect type of sampling rate. Type should be int or float.'),
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
