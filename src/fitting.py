from scipy import stats
import matplotlib.pyplot as plt

from nbinom_reparameterisation import nbinom_reparam, get_nbinom_params, method_of_moments, nbinom_reparam_cdf

def MLE(data, distribution, report=False):
    if distribution is nbinom_reparam:
        MOM_m, _ = method_of_moments(data, report)
        bounds = ((MOM_m, MOM_m),(0.01, 100)) # k smaller than 0.01 gives int errors, too large

        fitted_distribution = stats.fit(distribution, data, bounds)
        params = get_nbinom_params(fitted_distribution)
        
        if report:
            print("MAXIMUM LIKELIHOOD ESTIMATION mu={0:.04f}, k={1:.04f}".format(params.m, 
                                                                                 params.k))
        return fitted_distribution
    elif distribution is stats.geom:
        bounds = (0,1)
        fitted_distribution = stats.fit(distribution, data)
        if report:
            print("MAXIMUM LIKELIHOOD ESTIMATION p={0:.04f}".format(fitted_distribution.params.p))
        return fitted_distribution
    elif distribution is stats.poisson:
        # The parameter is just the sample mean... but let scipy do it
        bounds = [(0,12)]
        fitted_distribution = stats.fit(distribution, data, bounds)
        if report:
            print("MAXIMUM LIKELIHOOD ESTIMATION mu={0:.04f}".format(fitted_distribution.params.mu))
        return fitted_distribution

# See the notebook on calculating scaled AIC
def get_delta_AIC_c(data, dist, params):
    k = len(params)
    n = len(data)
    log_lik = dist.logpmf(data,*params).sum()
    AIC = 2*k - 2*log_lik
    AIC_c = AIC + (2*k*(k+1))/(n - k - 1)
    delta_AIC_c = AIC_c - min(data)
    return delta_AIC_c

def plot_goodness_of_fit(e_cdf, coaggregation_per_nymph, years_display, tick_species_display, host_species_display):
    # NBINOM DIST
    coaggregation_distribution_nbinom = MLE(coaggregation_per_nymph, nbinom_reparam)
    fitted_coaggregation_distribution_params = get_nbinom_params(coaggregation_distribution_nbinom)
    nbinom_cdf = nbinom_reparam_cdf(e_cdf.cdf.quantiles,
                                    m=fitted_coaggregation_distribution_params.m,
                                    k=fitted_coaggregation_distribution_params.k_reciprocal)
    plt.step(e_cdf.cdf.quantiles, nbinom_cdf, label="nbinom CDF")

    # GEOM DIST
    coaggregation_distribution_geom = MLE(coaggregation_per_nymph, stats.geom)
    geom_cdf = stats.geom.cdf(e_cdf.cdf.quantiles, p=coaggregation_distribution_geom.params.p)
    plt.step(e_cdf.cdf.quantiles, geom_cdf, label="geom CDF")
    
    # POISSON DIST
    coaggregation_distribution_poisson = MLE(coaggregation_per_nymph, stats.poisson)
    poisson_cdf = stats.poisson.cdf(e_cdf.cdf.quantiles, mu=coaggregation_distribution_poisson.params.mu)
    plt.step(e_cdf.cdf.quantiles, poisson_cdf, label="poisson CDF")
    
    # EMPERICAL CDF
    plt.step(e_cdf.cdf.quantiles, e_cdf.cdf.probabilities, label="emperical CDF")

    plt.legend()
    plt.title("eCDF and fitted distribution CDFs")
    #plt.title("Ticks: {} & Vertebrates: {}".format(tick_species_display, host_species_display))
    plt.grid(axis='x', color='0.95')
    plt.xlabel("x")
    plt.ylabel("CDF(x)")
    filename = "figs/CDF_compare_{}_{}_{}.png".format(years_display, tick_species_display, host_species_display).replace(" ", "").replace(",", "")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()
    
    # NBINOM DIST RESIDUALS
    nbinom_cdf = nbinom_reparam_cdf(e_cdf.cdf.quantiles, 
                                    m=fitted_coaggregation_distribution_params.m, 
                                    k=fitted_coaggregation_distribution_params.k_reciprocal)
    nbinom_cdf_residuals = abs(e_cdf.cdf.probabilities - nbinom_cdf)
    plt.plot(e_cdf.cdf.quantiles, nbinom_cdf_residuals, label="nbinom", marker=".", linewidth=0)
    
    # GEOM DIST RESIDUALS
    geom_cdf = stats.geom.cdf(e_cdf.cdf.quantiles, p=coaggregation_distribution_geom.params.p)
    geom_cdf_residuals = abs(e_cdf.cdf.probabilities - geom_cdf)
    plt.plot(e_cdf.cdf.quantiles, geom_cdf_residuals, label="geom", marker="s", linewidth=0)

    # POISSON DIST RESIDUALS
    poisson_cdf = stats.poisson.cdf(e_cdf.cdf.quantiles, mu=coaggregation_distribution_poisson.params.mu)
    poisson_cdf_residuals = abs(e_cdf.cdf.probabilities - poisson_cdf)
    plt.plot(e_cdf.cdf.quantiles, poisson_cdf_residuals, label="poisson", marker="p", linewidth=0)
    
    plt.legend()
    plt.title("Error: difference between eCDF and fitted CDFs")
    #plt.title("Tick: {} & Vertebrates: {}".format(tick_species_display, host_species_display))
    plt.grid(axis='x', color='0.95')
    plt.xlabel("x")
    plt.ylabel("|eCDF(x) - CDF(x)|")
    filename = "figs/CDF_errors_{}_{}_{}.png".format(years_display, tick_species_display, host_species_display).replace(" ", "").replace(",", "")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()

    poisson_scaled_AIC = get_delta_AIC_c(coaggregation_per_nymph, stats.poisson,        [coaggregation_distribution_poisson.params.mu])
    geom_scaled_AIC    = get_delta_AIC_c(coaggregation_per_nymph, stats.geom,           [coaggregation_distribution_geom.params.p])
    nbinom_scaled_AIC  = get_delta_AIC_c(coaggregation_per_nymph, nbinom_reparam,       [fitted_coaggregation_distribution_params.m,
                                                                                        fitted_coaggregation_distribution_params.k_reciprocal])

    print("AIC OF CO-AGGREGATION DISTRIBUTIONS")
    print("Poisson scaled AIC:", poisson_scaled_AIC)
    print("Geom scaled AIC:",    geom_scaled_AIC)
    print("NBinom scaled AIC:",  nbinom_scaled_AIC)

