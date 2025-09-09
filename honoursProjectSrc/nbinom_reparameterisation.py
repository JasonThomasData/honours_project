import sys
from scipy import stats
from math import factorial, gamma
from scipy.stats._distn_infrastructure import _ShapeInfo
import numpy as np

class NbinomReparameterised(stats.rv_discrete):

    def _pmf(self, x, m, k):
        
        pmf_values = []

        if isinstance(m, np.ndarray):
            m = m[0] # an array all the same for n, p
        if isinstance(k, np.ndarray):
            k = k[0]

        # Please see the Jupyter Notebook about the reparameterisation of the NBinom distribution for why using 1/k is the best choice
        k = 1/k
        
        # the standard method seems to use a np array, but I couldn't get it to work. So, verbosely use each element rather than a np array
        for x_i in x:
            if not isinstance(x, int):
                x_i = int(x_i)
            pmf_i = (1+m/k)**(-k) * (gamma(k+x_i)/(factorial(x_i)*gamma(k))) * (m/(m+k))**x_i
            pmf_values.append(pmf_i)
        return pmf_values
    
    def _argcheck(self, m, k):
        m_check = 0 < m
        k_check = 0 < k
        # there's no reason to assume that k would be less than one in this parameterisation
        return m_check & k_check

    def _param_info(self):
        # required to work with MLE in scipy.stats.fit
        return [_ShapeInfo("m", False, (0, np.inf), (True, False)),
                _ShapeInfo("k", False, (0, np.inf), (True, False)),
                _ShapeInfo("loc", False, (0, np.inf), (True, True))]
    
nbinom_reparam = NbinomReparameterised(name='NbinomReparameterised')

def nbinom_reparam_cdf(x, m, k):
    # I tried several times to add CDF methods to this distribution
    # So, let's just use this function outside the distribution's definition
    # Rather than derive the CDF formula for the reparameterised function, this is a list/vector of PMF values
    if isinstance(x, range):
        x = list(x)
    pmf_values = []
    last_x_i = -1 # this will make it so we always start the summation from 0
    for x_i in x:
        x_to_add = list(range(int(last_x_i+1), int(x_i+1)))
        pmf_value = np.sum(nbinom_reparam.pmf(x_to_add, m, k))
        pmf_values.append(pmf_value)
        last_x_i = x_i
    return np.cumsum(pmf_values)

def method_of_moments(data, report):
    mu = np.mean(data)
    sigma_squared = np.var(data) # Not the same sigma as used to calculate alpha
    m = mu
    k = mu**2/(sigma_squared - mu)
    if report:
        print("METHOD OF MOMENTS m={0:.04f}, k={1:.04f}".format(m, k))
        print("Mean: ", np.mean(data))
    return m, k

def get_nbinom_params(fitted_distribution):
    # Since we fit 1/k, we should now convert it back for when we need just k
    class Struct():
        def __init__(self):
            self.m = None
            self.k = None
            self.k_reciprocal = None
    params = Struct()
    params.m = fitted_distribution.params.m
    params.k = 1/fitted_distribution.params.k # use this where you need to use k itself
    params.k_reciprocal = fitted_distribution.params.k # use this where you need to use the PMF itself, since the PMF uses 1/k
    return params

