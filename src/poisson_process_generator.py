from scipy.stats import gamma, poisson

class PoissonProcessGenerator:
    # presets are to pregenerate this RV in tests for testable results
    def __init__(self, shape, scale, presets=[]):
        self.shape = shape
        self.scale = scale
        self.presets = presets

    def get(self):
        if len(self.presets) == 0:
            v = gamma.rvs(a=self.shape, scale=self.scale);
            return int(poisson.rvs(v,size=1)[0])
        else:
            return self.presets.pop(0)
