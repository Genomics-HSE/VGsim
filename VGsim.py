class Simulator():
	"""
	The class which creates and runs simulations.
	"""

	def __init__(number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, sampling_probability=False, memory_optimization=None):
		"""

		:param number_of_sites: the number of mutable sites
		:type number_of_sites: int

		:param populations_number: the number of populations (demes)
		:type populations_number: int

		:param number_of_susceptible_groups: the number of susceptible groups (groups with different immunity response)
		:type number_of_susceptible_groups: int

		:param seed: seed to generate simulation from. If **None**, then chosen at random
		:type seed: float or None

		:param sampling_probability: whether we set sampling probability as a share of recovered individuals (**True** value) or we will set it explicitly. Default is **False**.
		:type sampling_probability: bool or None

		:param memory_optimization: if True, then memory optimization is conducted (useful for large number of possible haplotypes)
		:type memory_optimization: bool or None
		"""
