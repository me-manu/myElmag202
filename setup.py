from numpy.distutils.core import Extension,setup

ext1 = Extension(name = 'elmagcore',
                 sources = ['./src/user202f2py.f90',
			     './src/aux202.f90',
			     './src/modules202.f90',
			     './src/init202.f90',
			     './src/elmag202.f90']
		 )

if __name__ == "__main__":
    setup(name = 'elmag',
          description       = "F2PY ELMAG202 Implementation",
          author            = "Manuel Meyer",
          author_email      = "manuel.meyer@fysik.su.se",
          version	    = "0.1",
          ext_modules = [ext1]
          )
# End of setup_example.py
