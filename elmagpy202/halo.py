"""
Some helper functions concerning the halo calculation.
"""

# --- Imports -------------------------------------------- #
import numpy as np
from scipy.interpolate import UnivariateSpline as USpline
import warnings
# -------------------------------------------------------- #

def calc_conf_radius(hist_etht, theta, conf, 
	start = 1,
	tmax = None, 
	t_lbounds = None, 
	emin = None, emax = None,
	e_lbounds = None,
	interp = True):
    """
    Calculate the containment radius of the halo component for all energies 
    for a given confidence level.

    Parameters
    ----------
    hist_etht:	(n,m,(k))-dim array, output of ELMAG simulation, 
		n: energies, m: angular sep., k: time delay. k dimension is optional.
    theta:	m-dim array, bin centers of angular seperation
    conf:	float, 0 < conf < 1, desired containment level

    kwargs
    ------
    emin:	float or None, minimum considered energy of halo photons in eV
    emax:	float or None, minimum considered energy of halo photons in eV
    e_lbounds:	(n+1)-dim array, left bin bounds of the energy in eV
    tmax:	float or None, maximum time delay of the cascade in years. 
    		If given, you need to provide the time delay bin bounds with 
		the t_lbounds keyword
    t_lbounds:	(j+1)-dim array, log10 of left bin bounds of the time delay
    interp:	bool, if true, use spline interpolation to find the containment 
    		radius. If false, use the nearest neighbour below the desired 
		containment level
    start:	int, 
    		first theta bin to be used (0: all photons are cascade photons)

    Returns
    -------
    float, containment radius in degrees
    """
    if conf < 0 or conf > 1:
	raise ValueError("Confidence level must be >= 0 and <= 1, not {0:.3f}".format(conf))

    if not tmax == None:
	idt	= np.where(t_lbounds <= np.log10(tmax)) # choose a maximum time delay
	if not len(idt):
	    raise ValueError("No cascade photons pass the chosen max. time cut of {0:.3e} years!".format(tmax))
	elif idt[0][-1] == 0:
	    raise ValueError("No cascade photons pass the chosen max. time cut of {0:.3e} years!".format(tmax))

	idt	= idt[0][-1] 
	hist	= hist_etht[:,start:,:idt].sum(axis = 2)	# index 1 in column 1: do not consider
							# primary photons, sum over all time delays
    elif len(hist_etht.shape) == 2:	# case of no time axis
	hist = hist_etht
    else:
	hist = hist_etht[:,start:,:].sum(axis = 2)

    if not emax == None:
	ide  = np.where(e_lbounds <= emax)
	if not len(ide):
	    raise ValueError("No cascade photons pass the chosen max. energy cut of {0:.3e} eV!".format(emax))
	elif ide[0][-1] == 0:
	    raise ValueError("No cascade photons pass the chosen max. energy cut of {0:.3e} eV!".format(emax))
	ide = ide[0][-1]
	hist = hist[ :ide if ide <= hist.shape[0] else None ,:]
    if not emin == None:
	ide  = np.where(e_lbounds >= emin)[0][0]
	hist = hist[ide:,:]

    if np.all(hist < 1e-10):	# no cascade photons
	warnings.warn("All entries smaller than 1e-10. No cascade photons? Retuning 0.",RuntimeWarning)
	return 0.
    hist = np.cumsum(hist.sum(axis = 0)) 
    hist /= hist[-1]	# CDF

    if interp:
	idmax = np.where(np.round(hist,5) == 1.)	# first entry where cdf is 1
	if not len(idmax) or not len(idmax[0]):
	    sh = slice(None)
	    st = slice(1,None)
	else:
	    idmax = idmax[0][0]
	    sh = slice(idmax+1 if idmax <= hist.shape[0] else None)
	    st = slice(start,idmax+2 if idmax+1 <= theta.shape[0] -1 else None)

	# interpolate inverted cdf
	# exclude first bin with primary emission
	try:
	    cdf_inv = USpline(hist[sh], np.log10(theta[st]), 
		    s = 0, k = 1, ext = 'raise')
	    if np.isnan(cdf_inv(conf)):
		# interpolation did not work, likely because there are not 
		# many cascade photons with the applied cuts. Using nearest bin instead.
		idx = np.where(hist <= conf)
		return (theta[start:])[idx[0][-1]]
	except:
		idx = np.where(hist <= conf)
		if not len(idx) or not len(idx[0]):
		    warnings.warn("No confidence level below {0:.3f}, minimum is {1:.3f}. Retruning 0".format(
			conf, hist[0]), RuntimeWarning)
		    return 0.
		return (theta[start:])[idx[0][-1]]
	return np.power(10.,cdf_inv(conf))

    else:
	idx = np.where(hist <= conf)
	if not len(idx):
	    warnings.warn("No confidence level below {0:.3f}, minimum is {1:.3f}. Retruning 0".format(
		conf, hist[0]), RuntimeWarning)
	    return 0.
	return (theta[start:])[idx[0][-1]]
