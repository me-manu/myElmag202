"""
Python wrapper for ELMAG code version 2.02
"""
__version__ = 0.1
__author__ = "Manuel Meyer // manuel.meyer@fysik.su.se"

# --- Imports --- #
import elmagcore
import numpy as np
import logging
import os
import time
from elmagpy202 import defaults
from os import path
# --------------- #
class Elmag(object):

    fname = __module__.replace('.','/') + '.py'

    @defaults.setDefault
    def __init__(self,model,**kwargs):
	"""
	Init the PyElmag class.

	Parameter
	---------
	model:		int, ebl model id,
			has to between 1 and 6

	kwargs
	------
	igmf:		float, egmf (in gauss), default: 1e-16
	cohlnth:	float, coherence length (in in Mpc), default: 1.
	gam1:		float, photon index 1, default: -1.7
	gam2:		float, photon index 2, default: -1.7
	nmax:		float, number of injected particles, default: 1e3
	ir_rescale:	float, ebl scale, default: 1
	th_jet:		float, jet opening angle in degree, default: 6
	a_smp:		float between 0 and 1, 
			sampling of the cascade, default: 0 (all particles are traced)
	ethr:		float, energy threshold for cascade simulation, eV, default: 1e8
	egmax:		float, maximal gamma energy, eV, default: 15e12
	egdelta:	float, gamma energy for delta distribution (at rest in Hubble flow), eV, default: 1e12
	emin:		float, low energy cutoff of injected photons, default 1e9
	ebreak:		float, break energy of injected spectrum, default: 15e12 (i.e. no break)
	log10t_min: 	float, minimum time in log10( dt/years )
	log10t_max: 	float, maximum time in log10( dt/years )
	log10theta_min:	float, minimum angular distance in log10( dtheta/deg ).
	log10theta_max:	float, maximum angular distance in log10( dtheta/deg )
	n_bine:		int, number of energy bins
	n_binth:	int, number of angular bins
	n_bint:		int, number of time bins
	tabledir:	string, directory with tabulated EBL files, default: "Tables/"
	"""

	self.model = model

	self.__dict__.update(kwargs)
	if self.model < 1 or self.model > 6:
	    raise ValueError("Uknown EBL model: model = {0:n}. Has to be between 1 and 6".format(
		self.model)
		)

	# save the current dir
	self.pwd = os.environ["PWD"]

	self.replaceZero = 1e-20 # float that replaces zeros in elmag output

	try:
	    self._e
	    logging.warning("Careful, you are trying to re-initialize ELMAG, this would crash python.")
	    self.set_params()
	    self.set_result_arrays()
	except AttributeError:
	    self._e = elmagcore # ugly work around?
	    self.set_tablefiles()
	    self.set_params()
	    self.set_result_arrays()
	    self._e.init(0,0)

	return

    def set_tablefiles(self):
	"""
	Set the names of the EBL density and redshift tables in a fortran readable way.
	"""
	tablefiles = [ path.join(self.tabledir,defaults.tables_n[self.model - 1]) ,
			path.join(self.tabledir,defaults.tables_z[self.model - 1]) ,
			path.join(self.tabledir,'redshift')
			]

	for it,t in enumerate(tablefiles):
	    if not it:
		fortran_table = self._e.user_variables.tablefile_n
	    elif it == 1:
		fortran_table = self._e.user_variables.tablefile_z
	    elif it == 2:
		fortran_table = self._e.user_variables.tablefile_redshift
	    # some checks
	    if len(t) > len(fortran_table):
		raise ValueError("String {0:s} has length {1:n} > {2:n}. Too long!".format(t,
						len(t), len(fortran_table))
			    )
	    if not path.isfile(t):
		raise ValueError("File {0:s} does not exist!".format(t))

	    # change paths
	    for i,s in enumerate(t):
		if not it:
		    self._e.user_variables.tablefile_n[i] = s
		elif it == 1:
		    self._e.user_variables.tablefile_z[i] = s
		elif it == 2:
		    self._e.user_variables.tablefile_redshift[i] = s

	    # fill up remaining characters with blanks:
	    for j in range(i+1,len(fortran_table)):
		if not it:
		    self._e.user_variables.tablefile_n[j] = ''
		elif it == 1:
		    self._e.user_variables.tablefile_z[j] = ''
		elif it == 2:
		    self._e.user_variables.tablefile_redshift[i] = s
	    logging.info('Set table file to {0:s}'.format(t))
	return

    def set_params(self,**kwargs):
	"""
	Update elmagcore parameters

	kwargs
	------
	Same as init function (you cannot change the EBL model).
	"""

	self.__dict__.update(kwargs)
	# ugly workaround:
	self._e.user_result.n_bine = self.n_bine
	self._e.user_result.n_bint = self.n_bint
	self._e.user_result.n_binth = self.n_binth

     	self._e.user_result.log10t_min = self.log10t_min
     	self._e.user_result.log10t_max = self.log10t_max
     	self._e.user_result.log10theta_min = self.log10theta_min
     	self._e.user_result.log10theta_max = self.log10theta_max

	self._e.user_variables.ethr = self.ethr
  	self._e.user_variables.egmax = self.egmax
  	self._e.user_variables.ir_rescale = self.ir_rescale
  	self._e.user_variables.cohlnth=3.086e21*1e3* self.cohlnth
  	self._e.user_variables.th_jet = self.th_jet
  	self._e.user_variables.a_smp = self.a_smp
	self._e.user_variables.igmf = self.igmf

  	self._e.user_variables.emin = self.emin
  	self._e.user_variables.ebreak = self.ebreak
  	self._e.user_variables.gam1 = self.gam1
  	self._e.user_variables.gam2 = self.gam2

  	self._e.user_variables.iseed = self.iseed

	# left bounds of theta and t bins in log10 representation
	# all events smaller than the zeroth bin bound end up in the first bin
	# all events larger than the -1-th bin bound end up in the last bin
	self.th_lbounds = np.linspace(self.log10theta_min,self.log10theta_max,self.n_binth + 1)
	self.t_lbounds = np.linspace(self.log10t_min,self.log10t_max,self.n_bint + 1)

	# bin centers
	self.th = np.power(10., 0.5 * (self.th_lbounds[1:] + self.th_lbounds[:-1]))
	self.t = np.power(10., 0.5 * (self.th_lbounds[1:] + self.th_lbounds[:-1]))
	
	# energy bounds and bin centers
	# in the original elmag output, elmag does not give the last bin of the histogram
	self.EeV_lbounds  = np.exp(np.linspace(np.log(self.ethr),np.log(self.egmax),self.n_bine))
	self.EeV = np.sqrt(self.EeV_lbounds[:-1] * self.EeV_lbounds[1:])	# bin centers

	return 

    def set_result_arrays(self):
	"""
	(Re-)set the result arrays to zeros. 
	"""
	self._e.user_result.hist_ge = np.zeros((self._e.user_result.n_bine,2))
	self._e.user_result.hist_eth = np.zeros((self._e.user_result.n_bine,
						self._e.user_result.n_binth))
	self._e.user_result.hist_et = np.zeros((self._e.user_result.n_bine,
						self._e.user_result.n_bint))
	self._e.user_result.hist_etht = np.zeros((self._e.user_result.n_bine,
						self._e.user_result.n_binth,
						self._e.user_result.n_bint))
	self.init_hist = np.zeros(self.EeV_lbounds.shape[0])
	return

    def initial_particle(self,distr= 'elmag',species = 0, egmax = 0., emin = 0.):
	"""
	Draw an initial particle energy. Sets e0, weight, and icq attributes.

	Returns
	-------
	tuple with energy, weight and particle species

	kwargs
	------
	distr:	str, determine from which distribution the particle is drawn.
		Options:
		    - 'elmag': use elmag's initial particle function (assumes broken power law)
		    - 'delta': delta distribution, energy will be equal to egdelta
	species: int, 0 for gamma ray, 1 for positron / electron
	egmax:	float, maximum possible energy of initial particle in eV. Default: equal to self.egmax
	emin:	float, minimum possible energy of initial particle in eV. Default: equal to self.emin
	"""
	egmax = self.egmax if egmax == 0. else egmax
	emin = self.emin if emin == 0. else emin
	if emin < self.emin:
	    raise ValueError("Requested emin={0:.3e} eV is below the emin set for ELMAG, {1:.3e} eV ".format(emin, self.emin))
	if egmax > self.egmax:
	    raise ValueError("Requested egmax={0:.3e} eV is above the egmax set for ELMAG, {1:.3e} eV ".format(egmax, self.egmax))

	if distr == 'elmag':
	    #self.e0, self.weight = self._e.initial_particle()
	    self.e0=emin*(egmax/emin)**self._e.psran() #energy: uniform in ln E 
	    if self.e0 < self.ebreak:
		#self.weight=(self.e0/self.ebreak)**(self.gam1+1.)*np.log(self.egmax/self.emin)
		self.weight=(self.e0/self.ebreak)**(self.gam1+1.)*np.log(egmax/emin)
	    else:
		#self.weight=(self.e0/self.ebreak)**(self.gam2+1.)*np.log(self.egmax/self.emin) 
		self.weight=(self.e0/self.ebreak)**(self.gam2+1.)*np.log(egmax/emin) 

	elif distr == 'delta':
	    self.e0, self.weight = self.egdelta, 1.
	else:
	    raise ValueError('{0:s} is unknown distribution'.format(distr))
	self.species = species
	return self.e0, self.weight, self.species

    def save_initial_particle(self,z):
	"""
	add initial particle to a histogram
	"""
	i = np.min( (self.n_bine, 
		int(np.log(self.e0/self.ethr)/np.log(self.egmax/self.ethr)*(self.n_bine-1)-np.log(1. + z)) + 1 )
		)
	i = np.max((i,1))
	self.init_hist[i-1] += self.weight * \
			      1./np.log(self.egmax/self.ethr)*(self.n_bine-1.)
			      #self.e0/np.log(self.egmax/self.ethr)*(self.n_bine-1.) / (1.+z)
	return

    def cascade(self,z):
	"""
	Run ELMAG to calculate the cascade.

	Parameters
	----------
	z:	float, source redshift
	"""
	self._e.cascade(self.species,self.e0,self.weight,z)
	return 

    def get_results(self):
	"""
	Collect the result histograms. Replace zeros with replaceZero attribute.
	"""
	# exclude the last energy bin, it's for energies > egmax
	self.hist_ge = self._e.user_result.hist_ge[:-1,:]
	self.hist_eth = self._e.user_result.hist_eth[:-1,:]
	self.hist_et = self._e.user_result.hist_et[:-1,:]
	self.hist_etht = self._e.user_result.hist_etht[:-1,:,:]
	self.init_hist = self.init_hist[:-1]
	
	self.hist_ge[self.hist_ge== 0.] = self.replaceZero * np.ones(np.sum(self.hist_ge == 0.))
	self.hist_eth[self.hist_eth== 0.] = self.replaceZero * np.ones(np.sum(self.hist_eth == 0.))
	self.hist_et[self.hist_et== 0.] = self.replaceZero * np.ones(np.sum(self.hist_et == 0.))
	self.hist_etht[self.hist_etht== 0.] = self.replaceZero * np.ones(np.sum(self.hist_etht == 0.))

	return

    def bin_by_bin_cascade(self, z, nmax, ibin = None, ebl_weight = 1):
	"""
	Run the cascade for each energy bin separately. 
	Save the results with the get_results function.

	Parameters
	----------
	z:	float, source redshift.
	nmax:	int, minimum number of injected particles per bin.

	kwargs
	------
	ibin:	int or None. If int, run only for this bin number where ibin = 0,...,self.nbin_e - 1
	ebl_weight: float. If > 1: increase bin sampling with optical depth tau. 
		    For each bin nmax * weight particles will be injected, where 
		    weight = floor(tau) + 1
		    weight is bounded by ebl_weight.
		    Requires eblstud package.
	"""
	if ebl_weight > 1:
	    from eblstud.ebl.tau_from_model import OptDepth as OD
	    tau = OD(model = defaults.eblmodel[self.model-1])

	    weights = np.floor(tau.opt_depth_array(z, self.EeV / 1e12)[0]) + 1
	    weights[weights > ebl_weight] = np.ones(np.sum(weights > ebl_weight)) * ebl_weight
	else:
	    weights = np.ones(self.EeV.shape)

	t = time.time()
	nmax_tot = 0

	for i,w in enumerate(weights):
	    if not ibin == None:
		if not i == ibin: continue
	    for j in range(int(w * nmax)):
		self.initial_particle(distr = 'elmag', emin=self.EeV_lbounds[i],
                          egmax = self.EeV_lbounds[i+1]) 
		self.weight /= w
		self.save_initial_particle(z)
		self.cascade(z)
		nmax_tot+=1

	logging.info("Done simulating {0:n} cascades. It took {1:.2f}s".format(nmax_tot,time.time() - t))

	self.get_results()

	return
