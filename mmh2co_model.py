from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import astropy.units as u
import scipy.ndimage as nd

#from pyspeckit.spectrum.units import SpectroscopicAxis 
t = Table.read('ph2cogrid.fits.gz',format='fits')

def H2COModel(t):
    Temps = t['Temperature'].data
    logN = np.log10(t['Column'].data)
    logn = np.log10(t['nH2'].data)
    logFWHM = np.log10(t['FWHM'].data)

    Tempbins = np.unique(Temps)
    logNbins = np.unique(logN)
    lognbins = np.unique(logn)
    logFWHMbins = np.unique(logFWHM)
    
    Tempaxis = interp.interp1d(Tempbins,np.arange(len(Tempbins)))
    logNaxis = interp.interp1d(logNbins,np.arange(len(logNbins)))
    lognaxis = interp.interp1d(lognbins,np.arange(len(lognbins)))
    logFHWMaxis = interp.interp1d(logFWHMbins,np.arange(len(logFWHMbins)))
        
    TempIdx = np.digitize(Temps,Tempbins,right=True)
    logNIdx = np.digitize(logN,logNbins,right=True)
    lognIdx = np.digitize(logn,lognbins,right=True)
    logFWHMIdx = np.digitize(logFWHM,logFWHMbins,right=True)

    
    Tex303 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    Tex303[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['Tex_303_202'].data

    Tex322 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    Tex322[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['Tex_322_221'].data

    Tex321 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    Tex321[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['Tex_321_220'].data

    tau303 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    tau303[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['tau_303_202'].data

    tau322 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    tau322[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['tau_322_221'].data

    tau321 = np.zeros(Tempbins.shape+logNbins.shape+\
                      lognbins.shape+logFWHMbins.shape)
    tau321[TempIdx,logNIdx,lognIdx,logFWHMIdx]=t['tau_321_220'].data

    def ModelFunction(Temperature = 31,logColumn = 13.5,
                      logDensity = 3.3, logwidth = 0.2):
        try:
            coords = np.array([[Tempaxis(Temperature),logNaxis(logColumn),
                            lognaxis(logDensity),logFHWMaxis(logwidth)]])
        except ValueError:
# If the value is out of bounds on an axis, ValueError
# is thrown so return some NaNs
            output = {'Tex303':np.nan,
                      'Tex322':np.nan,
                      'Tex321':np.nan,
                      'tau303':np.nan,
                      'tau322':np.nan,
                      'tau321':np.nan}
            return output
        
        output = {'Tex303':(nd.map_coordinates(Tex303,coords.T))[0],
                  'Tex322':(nd.map_coordinates(Tex322,coords.T))[0],
                  'Tex321':(nd.map_coordinates(Tex321,coords.T))[0],
                  'tau303':(nd.map_coordinates(tau303,coords.T))[0],
                  'tau322':(nd.map_coordinates(tau322,coords.T))[0],
                  'tau321':(nd.map_coordinates(tau321,coords.T))[0]}
        return output

    return ModelFunction
