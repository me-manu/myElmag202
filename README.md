#myElmag202

A python implementation of the ELMAG 202 Monte-Carlo code. ELMAG calculates the electromagnetic cascade initiated by high energy gammarays. The official code website can be found [here](http://elmag.sourceforge.net/). The reference paper of the code can be found [here](http://adsabs.harvard.edu/abs/2012CoPhC.183.1036K).

This implementation provides all necessary files to install and run ELMAG from python. It relies on the f2py fortran to python conversion included in numpy. The implementation here also provides an additional model for the extragalactic background light (EBL), namely the one from [Dominguez et al. (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.410.2556D).

##Prerequisites

You will need numpy with f2py. 
Optionally, you can also install the package [eblstud](https://github.com/me-manu/eblstud) for the interpolation of EBL models. However, you only need it if you want to perform bin-by-bin simulations with ELMAG where each energy bin is weighted with the optical depth (in the `bin_by_bin_cascade` function in elmag202.py).

##Installation

After downloading the code, run 
```
python setup.py build
python setup.py install
```
The first line will compile the ELMAG source code and put the output in the `./build/` directory in current working directory. The second line will install the python egg. You might need root privileges to run this line. Alternatively, you can install the egg in a custom location by adding the flag `--prefix=/directory/with/writing/priviliges/`. In this case, make sure to add this directory to your python path. Assuming you installed ELMAG with python 2.7, type this (in bash, or add to your .bashrc file):

```
export PYTHONPATH=/prefix/dir/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=/directory/with/myELMAG202/:$PYTHONPATH
```

The installation has been successfully tested with:
- Python 2.7.5 [GCC 4.1.2 20080704 (Red Hat 4.1.2-54)] on linux2; gfortran gcc version 4.4.7 
- Python 2.7.10 [GCC 4.2.1 Compatible Apple LLVM 7.0.0 (clang-700.0.59.5)] on darwin; gcc version 5.3.0

##Testing

You can test the implementation by running the included ipython notebook.

