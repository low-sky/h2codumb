from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import astropy.units as u
from pyspeckit.spectrum.units import SpectroscopicAxis 
t = Table.read('ph2cogrid.fits.gz',format='fits')

execfile('h2co_mm.py')

Temps = t['Temperature']
logN = np.log10(t['Column'])
logn = np.log10(t['nH2'])
FWHM = t['FWHM']
axes =  np.array([logN.ravel(),logn.ravel(),Temps.ravel(),FWHM.ravel()]).T

Tex303 = interp.LinearNDInterpolator(axes,t['Tex_303_202'])
Tex322 = interp.LinearNDInterpolator(axes,t['Tex_322_221'])
Tex321 = interp.LinearNDInterpolator(axes,t['Tex_321_220'])
tau303 = interp.LinearNDInterpolator(axes,t['tau_303_202'])
tau322 = interp.LinearNDInterpolator(axes,t['tau_322_221'])
tau321 = interp.LinearNDInterpolator(axes,t['tau_321_220'])

bundle = (Tex303,Tex322,Tex321,tau303,tau322,tau321)

nu = SpectroscopicAxis(np.linspace(218e9,219e9,1000)*u.Hz)
spec = h2co_mm_radex(nu,gridbundle=bundle,verbose=True,Temperature=60)
plt.plot(nu.as_unit('GHz'),spec)
plt.show()
