from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import astropy.units as u
from pyspeckit.spectrum.units import SpectroscopicAxis 
from pyspeckit.spectrum.models.h2co_mm import h2co_mm_radex

import mmh2co_model as model

t = Table.read('ph2cogrid.fits.gz',format='fits')

bundle = model.H2COModel(t)

nu = SpectroscopicAxis(np.linspace(218e9,219e9,1000)*u.Hz)
spec = h2co_mm_radex(nu,gridbundle=bundle,verbose=True,Temperature=15,width=1.0)
plt.plot(nu.as_unit('GHz'),spec)
plt.show()
