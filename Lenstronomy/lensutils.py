import numpy as np
import math
import emcee
import sys
import pickle
from lenstronomy.Cosmo.kde_likelihood import KDELikelihood
from astropy.cosmology import FlatLambdaCDM, FlatwCDM, LambdaCDM


class StrongLensSystem(object):
    """
	This is a parent class, common to all lens modeling code outputs.
	It stores the "physical" parameters of the lens (name, redshifts, ...)
	"""

    def __init__(self, name, zlens, zsource, longname=None):
        self.name = name
        self.zlens = zlens
        self.zsource = zsource
        self.longname = longname

    def __str__(self):
        return "%s\n\tzl = %f\n\tzs = %f" % (self.name, self.zlens, self.zsource)


class GLEELens(StrongLensSystem):
    """
	This class takes the output of GLEE (Ddt distribution) from which it
    evaluates the likelihood of a Ddt (in Mpc) predicted in a given cosmology.

	The default likelihood follows a skewed log-normal distribution.
    No other likelihoods have been implemented so far.
	"""

    def __init__(
        self,
        name,
        zlens,
        zsource,
        loglikelihood_type="sklogn_analytical",
        mu=None,
        sigma=None,
        lam=None,
        explim=100.0,
        longname=None,
    ):

        StrongLensSystem.__init__(
            self, name=name, zlens=zlens, zsource=zsource, longname=longname
        )
        self.mu = mu
        self.sigma = sigma
        self.lam = lam
        self.explim = explim
        self.loglikelihood_type = loglikelihood_type
        self.init_loglikelihood()

    def sklogn_analytical_likelihood(self, ddt):
        """
		Evaluates the likelihood of a time-delay distance ddt (in Mpc) against
        the model predictions, using a skewed log-normal distribution.
		"""
        # skewed lognormal distribution with boundaries
        if (ddt < self.lam) or (
            (-self.mu + math.log(ddt - self.lam)) ** 2 / (2.0 * self.sigma ** 2)
            > self.explim
        ):
            return -np.inf
        else:
            llh = math.exp(
                -((-self.mu + math.log(ddt - self.lam)) ** 2 / (2.0 * self.sigma ** 2))
            ) / (math.sqrt(2 * math.pi) * (ddt - self.lam) * self.sigma)

            if np.isnan(llh):
                return -np.inf
            else:
                return np.log(llh)

    def init_loglikelihood(self):
        if self.loglikelihood_type == "sklogn_analytical":
            self.loglikelihood = self.sklogn_analytical_likelihood
        else:
            assert ValueError("unknown keyword: %s" % loglikelihood_type)
            # if you want to implement other likelihood estimators, do it here


class LenstronomyLens(StrongLensSystem):
    """
	This class takes the output of Lenstronomy (Dd versus Ddt distributions)
    from which it evaluates the likelihood of a Dd versus Ddt (in Mpc)
    predicted in a given cosmology.

	The default likelihood follows the KDE log-normal distribution implemented
    in Lenstronomy. You can change the type of kernel used. No other likelihoods
    have been implemented so far.
	"""

    def __init__(
        self,
        name,
        zlens,
        zsource,
        ddt_vs_dd_samples,
        longname=None,
        loglikelihood_type="kde",
        kde_type="scipy_gaussian",
    ):
        StrongLensSystem.__init__(
            self, name=name, zlens=zlens, zsource=zsource, longname=longname
        )

        self.ddt_vs_dd = ddt_vs_dd_samples
        self.loglikelihood_type = loglikelihood_type
        self.kde_type = kde_type
        self.init_loglikelihood()

    def kdelikelihood(self, kde_type, bandwidth=2):
        """
		Evaluates the likelihood of a angular diameter distance to the deflector
        Dd (in Mpc) versus its time-delay distance Ddt (in Mpc) against the
        model predictions, using a loglikelihood sampled
        from a Kernel Density Estimator.
		"""
        self.ddt = self.ddt_vs_dd["ddt"]
        self.dd = self.ddt_vs_dd["dd"]
        KDEl = KDELikelihood(
            self.dd.values, self.ddt.values, kde_type=kde_type, bandwidth=bandwidth
        )
        return KDEl.logLikelihood

    def init_loglikelihood(self):
        if self.loglikelihood_type == "kde":
            self.loglikelihood = self.kdelikelihood(kde_type=self.kde_type)
        else:
            assert ValueError("unknown keyword: %s" % self.loglikelihood_type)
            # if you want to implement other likelihood estimators, do it here


def log_prior(theta, cosmology):
    """
	Return flat priors on the cosmological parameters - hardcoded boundaries.

    Parameters:
    -----------
	theta:
        list of floats, folded cosmological parameters.
	cosmology:
        string, keyword indicating the choice of cosmology to work with.
	"""
    if cosmology == "FLCDM":
        h0, om = theta
        if 0.0 <= h0 <= 150.0 and 0.05 <= om <= 0.5:
            return 0.0
        else:
            return -np.inf

    elif cosmology == "FwCDM":
        h0, om, w = theta
        if 0.0 <= h0 <= 150.0 and 0.05 <= om <= 0.5 and -2.5 <= w <= 0.5:
            return 0.0
        else:
            return -np.inf

    elif cosmology == "oLCDM":
        h0, om, ok = theta
        if 0.0 <= h0 <= 150.0 and 0.05 <= om <= 0.5 and -0.5 <= ok <= 0.5:
            # make sure that Omega_DE is not negative...
            if 1.0 - om - ok <= 0:
                return -np.inf
            else:
                return 0.0
        else:
            return -np.inf


def log_like_add(lens, cosmo):
    """
	Computes the relevant angular diameter distance(s) of a given lens in a
    given cosmology, and evaluate its/their joint likelihood against
    the same modeled distances of the lens.

    Parameters:
    -----------
	lens:
        either a GLEELens or LenstronomyLens instance.
	cosmo:
        an astropy cosmology object. 
	"""
    dd = cosmo.angular_diameter_distance(z=lens.zlens).value
    ds = cosmo.angular_diameter_distance(z=lens.zsource).value
    dds = cosmo.angular_diameter_distance_z1z2(z1=lens.zlens, z2=lens.zsource).value
    ddt = (1.0 + lens.zlens) * dd * ds / dds

    if isinstance(lens, GLEELens):
        return lens.loglikelihood(ddt)

    elif isinstance(lens, LenstronomyLens):
        return lens.loglikelihood(dd, ddt)

    else:
        sys.exit("I don't know what to do with %s, unknown instance" % lens)


def log_prob_ddt(theta, lenses, cosmology):
    """
	Compute the likelihood of the given cosmological parameters against the
	modeled angular diameter distances of the lenses.

    Parameters
    ----------
	theta: list
        loat folded cosmological parameters.
	lenses: list
        lens objects (currently either GLEELens or LenstronomyLens).
	cosmology: string
        keyword indicating the choice of cosmology to work with.
	"""

    lp = log_prior(theta, cosmology)
    if not np.isfinite(lp):
        return -np.inf
    else:
        logprob = lp
        if cosmology == "FLCDM":
            h0, om = theta
            cosmo = FlatLambdaCDM(H0=h0, Om0=om)
        elif cosmology == "FwCDM":
            h0, om, w = theta
            cosmo = FlatwCDM(H0=h0, Om0=om, w0=w)
        elif cosmology == "oLCDM":
            h0, om, ok = theta
            # assert we are not in a crazy cosmological situation that prevents
            # computing the angular distance integral
            if np.any(
                [
                    ok * (1.0 + lens.zsource) ** 2
                    + om * (1.0 + lens.zsource) ** 3
                    + (1.0 - om - ok)
                    <= 0
                    for lens in lenses
                ]
            ):
                return -np.inf
            else:
                cosmo = LambdaCDM(H0=h0, Om0=om, Ode0=1.0 - om - ok)
        else:
            raise ValueError("I don't know the cosmology %s" % cosmology)
        
        for lens in lenses:
            logprob += log_like_add(lens=lens, cosmo=cosmo)
        
        return logprob


def sample_params(
    lenses, cosmology, nwalkers=32, nsamples=20000, save=True, filepath="temp.pkl"
):
    """
	High-level wrapper around the log_prob_ddt() function.
    Explore the cosmological parameters space and return their likelihood
    evaluated against the modeled angular diameter distances
    of (multiple) lens system(s).

    Parameters
    ----------
	lenses: np.list
        lens objects (currently either GLEELens or LenstronomyLens).
	cosmology: np.str
        keyword indicating the choice of cosmology to work with.
	nwalkers: np.int
        number of emcee walkers used to sample the parameters space.
	nsamples: np.int
        number of samples for an MCMC chain to converge. 
	    Make sure these are larger than the autocorrelation time!
	save: np.boolean
        if True the combined, flattened chain is saved in filepath
	filepath: np.string
        path of where the output chain is saved.
	"""

    # Our starting point is a "decent" solution. You might want to check it does not impact the results, but unless you start from really crazy values, it shouldn't. We slightly randomize the starting point for each walker.
    if cosmology == "FLCDM":
        startpos = [72, 0.3]  # H0, Om
    elif cosmology == "FwCDM":
        startpos = [72, 0.3, -1]  # H0, Om, w
    elif cosmology == "oLCDM":
        startpos = [72, 0.3, 0.0]  # H0, Om, Ok
    else:
        raise ValueError("I don't know the cosmology %s" % cosmology)

    pos = startpos + 1e-4*np.random.randn(nwalkers, len(startpos))
    nwalkers, ndim = pos.shape

    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, log_prob_ddt, args=([lenses, cosmology])
    )
    sampler.run_mcmc(pos, nsamples, progress=True)

    tau = sampler.get_autocorr_time()
    print("Autocorrelation time: ", tau)
    flat_samples = sampler.get_chain(discard=1000, thin=15, flat=True)

    if save:
        pkl_file = open(filepath, "wb")
        pickle.dump(flat_samples, pkl_file, protocol=-1)
        pkl_file.close()

    return flat_samples


def readpickle(filepath):
    """
	Small utility function to read pickle files written by the sample_params
    function above.

	param filepath: path of the pickle file to load.
	"""
    pkl_file = open(filepath, "rb")
    try:
        obj = pickle.load(pkl_file, encoding="bytes")
    except:
        obj = pickle.load(pkl_file)
    pkl_file.close()
    return obj
