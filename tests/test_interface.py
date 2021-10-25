from numpy.testing import assert_allclose

from VGsim import Simulator

def test_interface():
    assert_allclose(2, 3, atol=1e-14)     # this fails, obviously
