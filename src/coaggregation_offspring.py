import numpy as np

from nbinom_reparameterisation import nbinom_reparam

# This is the bernoulli trial described in thesis
def infected_larva_becomes_infectious_nymph(alpha):
    return bool(np.random.uniform(low=0.0, high=1.0, size=1)[0] <= alpha)

# Described as I_i in thesis (binomial experiment)
def number_of_offspring_one_nymph(alpha, larvae_count):
    offspring = 0
    for _ in range(1,larvae_count+1):
        if infected_larva_becomes_infectious_nymph(alpha):
            offspring += 1
    return offspring

# Described as L_i in thesis
def get_coaggregation_per_nymph(dataframe, ricinus, trianguliceps):
    coaggregation_per_nymph = np.array([], dtype=np.int64)
    for i, row in dataframe.iterrows():
        nymph_count = 0
        larvae_count = 0
        if ricinus:
            nymph_count += row["ricinusN"]
            larvae_count += row["ricinusL"]
        if trianguliceps:
            nymph_count += row["triangulicepsN"]
            larvae_count += row["triangulicepsL"]
        larvae_coaggregating = np.repeat(larvae_count, nymph_count)
        coaggregation_per_nymph = np.append(coaggregation_per_nymph, larvae_coaggregating)
    return coaggregation_per_nymph

def get_aggregation_per_vertebrate(dataframe, ricinus, trianguliceps):
    aggregation_per_verterbrate = []
    for i, row in dataframe.iterrows():
        nymph_count = 0
        larvae_count = 0
        if ricinus:
            nymph_count += row["ricinusN"]
            larvae_count += row["ricinusL"]
        if trianguliceps:
            nymph_count += row["triangulicepsN"]
            larvae_count += row["triangulicepsL"]
        aggregation_per_verterbrate.append(int(nymph_count+larvae_count))
    return aggregation_per_verterbrate

def get_coaggregation_per_vertebrate_perkins(dataframe, ricinus, trianguliceps):
    # For the Perkins calculation, I could find no information about what a co-feeding group is, here's my best attempt
    coaggregation_per_verterbrate = []
    for i, row in dataframe.iterrows():
        nymph_count = 0
        larvae_count = 0
        if ricinus:
            nymph_count += row["ricinusN"]
            larvae_count += row["ricinusL"]
        if trianguliceps:
            nymph_count += row["triangulicepsN"]
            larvae_count += row["triangulicepsL"]

        if nymph_count > 0 and larvae_count > 0:
            coaggregation_per_verterbrate.append(int(nymph_count+larvae_count))
        else:
            coaggregation_per_verterbrate.append(int(0))
    return coaggregation_per_verterbrate

def get_coaggregation_per_vertebrate_JR(dataframe, ricinus, trianguliceps):
    # In thesis, b_n*b_l
    coaggregation_per_verterbrate = []
    for i, row in dataframe.iterrows():
        nymph_count = 0
        larvae_count = 0
        if ricinus:
            nymph_count += row["ricinusN"]
            larvae_count += row["ricinusL"]
        if trianguliceps:
            nymph_count += row["triangulicepsN"]
            larvae_count += row["triangulicepsL"]
        coaggregation_per_verterbrate.append(nymph_count*larvae_count)
    return coaggregation_per_verterbrate

def get_larval_burden_per_vertebrate(dataframe, ricinus, trianguliceps):
    larval_burden_per_vertebrate = []
    for i, row in dataframe.iterrows():
        larvae_count = 0
        if ricinus:
            larvae_count += row["ricinusL"]
        if trianguliceps:
            larvae_count += row["triangulicepsL"]

        larval_burden_per_vertebrate.append(larvae_count)
    return larval_burden_per_vertebrate

def simulate_offspring_data(alpha, fitted_coaggregation_distribution_params, simulation_size, report=True):
    if report:
        print("alpha = {0:.04f}".format(alpha))
    offspring_simulated = np.array([], dtype=np.int64)
    random_larvae_counts = nbinom_reparam.rvs(m=fitted_coaggregation_distribution_params.m,
                                              k=fitted_coaggregation_distribution_params.k_reciprocal,
                                              size=simulation_size)
    for larvae_count in random_larvae_counts:
        offspring_for_multiple_nymphs = number_of_offspring_one_nymph(alpha, larvae_count)
        offspring_simulated = np.append(offspring_simulated, offspring_for_multiple_nymphs)

    return offspring_simulated

