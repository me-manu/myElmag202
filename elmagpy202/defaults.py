# --- Imports --- #
import logging
import sys
# --------------- #

def setDefault(func):
    """
    Read in default keywords of the simulation and pass to function
    """
    def init(*args, **kwargs):
	for k in defaults.keys():
	    kwargs.setdefault(k,defaults[k])
	return func(*args, **kwargs)
    return init

defaults = {
    "log10t_min":     -3.00,
    "log10t_max":     9.00,
    "log10theta_min": -8.00,
    "log10theta_max": 1.30,
    "n_bine":         81,
    "n_binth":        93,
    "n_bint":         48,
    "model":          6,
    "ethr":           1e+8,
    "egmax":          15e+12,
    "egdelta":        1e+12,
    "ir_rescale":     1.00, 
    "cohlnth":        1.00, 
    "th_jet":         6.00, 
    "a_smp":          0.00,
    "igmf":           1e-16,
    "emin":           1.e+9,
    "ebreak":         1.5e+13,
    "gam1":           -1.70, 
    "gam2":           -1.70,
    "iseed":           1,
    "tabledir":	      "Tables/"
}

tables_n = ['n_bestfit10.dat',
	'n_lowerlimit10.dat',
	'n_Fra.dat',
	'n_Finke.dat',
	'n_Gil.dat',
	'n_Dom.dat'
]

tables_z = ['z-IR.dat',
    'z-IR.dat',
    'z-IR.dat',
    'z-IR.dat',
    'z-IR_Gil.dat',
    'z-IR_Dom.dat'
]

eblmodel = ['kneiske',
	'kneiske',
	'franceschini',
	'finke',
	'gilmore',
	'dominguez'
]

def init_logging(level, color = False):
    """
    Setup logger.

    Parameters 
    ----------
    level:	string, level of logging: DEBUG,INFO,WARNING,ERROR. (default: INFO).

    kwargs
    ------
    color:	bool, if true, enable colored output for bash output

    Notes
    -----
    for color see
	stackoverflow.com/questions/384076/how-can-i-color-python-logging-output
	https://wiki.archlinux.org/index.php/Color_Bash_Prompt
    """
    for handler in logging.root.handlers[:]:
	logging.root.removeHandler(handler)
    if level.upper() == 'INFO':
	level = logging.INFO
    elif level.upper() == 'DEBUG':
	level = logging.DEBUG
    elif level.upper() == 'WARNING':
	level = logging.WARNING
    elif level.upper() == 'ERROR':
	level = logging.ERROR


    if color:
	logging.basicConfig(level=level,stream = sys.stderr, format='\033[0;36m%(filename)10s:\033[0;35m%(lineno)4s\033[0;0m --- %(levelname)7s: %(message)s')
	logging.addLevelName( logging.DEBUG, "\033[1;32m%s\033[1;0m" % logging.getLevelName(logging.DEBUG))
	logging.addLevelName( logging.INFO, "\033[1;36m%s\033[1;0m" % logging.getLevelName(logging.INFO))
	logging.addLevelName( logging.WARNING, "\033[1;31m%s\033[1;0m" % logging.getLevelName(logging.WARNING))
	logging.addLevelName( logging.ERROR, "\033[1;41m%s\033[1;0m" % logging.getLevelName(logging.ERROR))
    else:
	logging.basicConfig(level=level,stream = sys.stderr, format='%(filename)10s:%(lineno)4s --- %(levelname)7s: %(message)s')

    return
