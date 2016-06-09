"""
Funtions to initiate and manipulate the bin-by-bin calculations with elmag.
"""

# --- Imports --------- #
import numpy as np
from scipy.integrate import simps
from elmagpy202 import halo
# --------------------- #

def calc_spectral_weights(EeV_centers, EeV_max, index_old, index_new):
    """
    Calculate weights for the cascade spectrum for each energy bin 
    to change the injection spectrum to an arbitrary spectrum.
    Only works for bin-by-bin simulations with alpha_smp = 0

    Parameters
    ----------
    EeV_centers:	n-dim array with energy centers of the injected spectrum.
			In units of eV.
    EeV_max:		float,
    			maximum energy used for simulation (egbreak in simulation).
			In units of eV.
    index_old:		float,
    			spectral index of simulated spectrum.
    index_new:		float,
    			desired spectral index. 

    Returns
    -------
    n-dim array with with weights. 
    """
    # set up energy integration array (in true energies)
    reweight = np.power(EeV_centers / EeV_max, index_new - index_old)
    return reweight

def apply_time_cut(bbbObs, t_lbounds, tmax):
    """
    Discard cascade photons that arrive later than a maximal time/

    Parameters
    ----------
    bbbObs:	(n , m , i , j) dimensional array with observed fluxes. 
		Axis are: true energy, obs. energy, angular separation, time delay
    t_lbounds:	j + 1 dim array with left bin boundaries of time delay in log10(years)
    tmax	float, maximum time delay in years

    Returns
    -------
    n x m x i dim array. All j bins with t > maxt are removed, and subseqeuntly it is summed over j axis
    """
    # find time bin for which left bin bound is smaller than tmax.
    idt = np.where(t_lbounds <= np.log10(tmax))
    if not len(idt):
	raise ValueError("No cascade photons pass the chosen max. time cut of {0:.3e} years!".format(tmax))
    elif idt[0][-1] == 0:
	raise ValueError("No cascade photons pass the chosen max. time cut of {0:.3e} years!".format(tmax))

    idt	= idt[0][-1] 
    return bbbObs[:,:,:,:idt].sum(axis = 3)

def calc_cont_radius_2d(bbbObsTimeCut,e_lbounds, theta, conf = 0.68, interp = False):
    """
    Calculate the containment radius for a given confidence level for each observed energy bin
    and each true energy bin. True energy bins are treated as maximum injection energies.

    Parameters
    ----------
    bbbObsTimeCut:	(n,m,k)-dim array with observed fluxes with a max time cut applied.
    			n: true energy, m: observed energy, k: angular separation. n should equal m.
    e_lbounds:		(n+1) = (m+1) dimensional array with left ebergy bin boundaries in eV
    theta:	k-dim array, bin centers of angular seperation in degrees

    kwargs
    ------
    conf:	float, 0 <= conf <= 1, confidence level of containment radius.
    		default: 0.68
    interp:	bool, if true, use spline interpolation to find the containment 
    		radius. If false, use the nearest neighbour below the desired 
		containment level

    Returns
    -------
    (n,m)-dim array with containment radii in degrees.
    """
    r = np.empty(bbbObsTimeCut.shape[:2])
    for i,et in enumerate(e_lbounds[1:]):	# loop over true energies
	# sum fluxes over true energies. This increases the maximum injection energy.
	# do not consider the primary gamma ray energies
	efluxSum = bbbObsTimeCut[:i+1 if i+1 < bbbObsTimeCut.shape[0] else None,:,1:].sum(axis = 0)
	for j,eo in enumerate(e_lbounds[1:]): # loop over observed energies
	    # calculate the containment radius for each observed energy bin 
	    r[i,j] = halo.calc_conf_radius(efluxSum, theta, conf, 
	    	tmax = None,
		emin = e_lbounds[j], 
		emax = eo,
		e_lbounds = e_lbounds,
		interp = interp)
	    #print "{0:.3e}   {1:.3e}   {2:.3f}".format(et / 1e12,eo/1e9,r[i,j])
    return r


#def calc_conf_radius():
