#!/usr/bin/env python

"""

from scipy import stats
from scipy.optimize import differential_evolution
import numpy as np

data = [0,0,0,1,2,4,11]

shape = 1.5
bounds = [(shape, shape),
          (0, np.inf)]

def optimizer(fun, bounds, *, integrality):
    print("opt function used")
    rng = np.random.default_rng
    return differential_evolution(fun, bounds, strategy='best2bin',
                                  rng=rng, integrality=[False,False,False])
								  
nbinom_fit = stats.fit(stats.nbinom, data, bounds, optimizer=optimizer)

print("MAXIMUM LIKELIHOOD ESTIMATION n={0:.04f}, p={1:.04f}".format(nbinom_fit.params.n, 
																	nbinom_fit.params.p))
                                                                    
"""

"""
from scipy import stats
from scipy.optimize import differential_evolution
import numpy as np

rng = np.random.default_rng()

def optimizer(fun, bounds, *, integrality):
    return differential_evolution(fun, bounds, strategy='best2bin',
                                  rng=rng, integrality=[False, False, False])

data = [0,0,0,1,2,4,11]
bounds = [(1.5, 1.55), (0, 1)]
res4 = stats.fit(stats.nbinom, data, bounds, optimizer=optimizer)
res4.params
"""

from scipy import stats
from scipy.optimize import differential_evolution
import numpy as np
from scipy.stats._discrete_distns import nbinom_gen, _ShapeInfo


class new_nbinom(nbinom_gen):
    def _shape_info(self):
        # in Shapeinfo we specify 'n' as non-integer
        return [_ShapeInfo("n", False, (0, np.inf), (True, False)),
                _ShapeInfo("p", False, (0, 1), (True, True))]
        

def optimizer(fun, bounds, *, integrality):
    return differential_evolution(fun, 
                                  bounds, 
                                  strategy='best2bin',
                                  integrality=[False, False, False])

data = [0,0,0,1,2,4,11]
bounds = [(0.1, 1.55), (0, 1)]
distr = new_nbinom(name = 'nbinom')

res4 = stats.fit(distr, data, bounds, optimizer=optimizer)
print(res4)
