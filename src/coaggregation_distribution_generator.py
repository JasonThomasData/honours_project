import numpy as np

class CoaggregationDistributionGenerator:
    def __init__(self, c, v, sigma, coaggregation_distribution):
		# v: the probability of transmission between two ticks
		# c: the contact rate, or, the probability that a nymph and larvae co-aggregate close enough for transmission to occur
		# sigma: the probability that an infected larvae survives to become an infections nymph
		self.offspring_threshold = v*c*sigma
        self.coaggregation_distribution = coaggregation_distribution

	def is_offspring(self):
		return bool(np.random.uniform(low=0.0, high=1.0, size=1) < self.offspring_threshold)

    def get(self):
		offspring = 0
		larvae_coaggregated_with_one_nymph = random.choice(self.coaggregation_distribution, size=1, replace=True)
		for _ in range(larvae_coaggregated_with_one_nymph):
			if self.is_offspring():
				offspring += 1
		return offspring
