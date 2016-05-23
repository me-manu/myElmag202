"""
Some helper functions concerning the halo calculation.
"""

# --- Imports -------------------------------------------- #
import numpy as np
from scipy.interpolate import UnivariateSpline as USpline
# -------------------------------------------------------- #

def calc_conf_radius(hist_etht, theta, conf, 
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
    hist_etht:	(n,m,k)-dim array, output of ELMAG simulation, 
		n: energies, m: angular sep., k: time delay
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

    Returns
    -------
    float, containment radius in degrees
    """
    if not tmax == None:
	idt	= np.where(t_lbounds <= np.log10(tmax))[0][-1] # choose a maximum time delay
	hist	= hist_etht[:,1:,:idt+1].sum(axis = 2)	# index 1 in column 1: do not consider
							# primary photons, sum over all time delays
    else:
	hist = hist_etht[:,1:,:].sum(axis = 2)

    if not emax == None:
	ide  = np.where(e_lbounds <= emax)[0][-1] # choose a maximum time delay
	hist = hist[ :ide+1 if ide <= hist.shape[0] else None ,:]
    if not emin == None:
	ide  = np.where(e_lbounds >= emin)[0][0] # choose a maximum time delay
	hist = hist[ide:,:]

    hist = np.cumsum(hist.sum(axis = 0)) / hist.sum()

    if interp:
	idmax = np.where(np.round(hist,5) == 1.)[0][0]

	cdf_inv = USpline(hist[:idmax+1 if idmax <= hist.shape[0] else None ], 
		    np.log10(theta[1:idmax+2 if idmax+1 <= theta.shape[0] -1 else None ]), 
		    s = 0, k = 1, ext = 'raise')
	return np.power(10.,cdf_inv(conf))

    else:
	idx = np.where(hist <= conf)[0][-1]
	return (theta[1:])[idx]
