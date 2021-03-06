from astropy.table import Table
import pyradex
import numpy as np

R = pyradex.Radex(column=1e16,abundance=1e-4,species='ph2co-h2')

fortho = 0.75

nFWHM = 5
FHWMmin = 0.5
FHWMmx = 5

nDens = 21
nlower = 2
nupper = 6

nCol = 41
Nlower = 10
Nupper = 16

nTemp = 31
Tlower = 10
Tupper = 300

Temps = np.logspace(1,2.5,nTemp)
Cols = 1e1**np.linspace(Nlower,Nupper,nCol)
Densities = 1e1**(np.linspace(nlower,nupper,nDens))
FWHM = np.logspace(np.log10(0.5),np.log10(5),nFWHM)


outtable = Table(names = ['Tex_303_202','Tex_322_221','Tex_321_220',
                          'tau_303_202','tau_322_221','tau_321_220',
                          'Temperature','Column','nH2','FWHM'])

for T in Temps:
    for N in Cols:
        for n in Densities:
            for dV in FWHM:
                Tlvg = R(collider_densities={'oH2':n*fortho,'pH2':(1-fortho)*n}, column=N, abundance = 1e-9, species='ph2co-h2',temperature=T,deltav=dV)
                outtable.add_row()
                outtable[-1]['Tex_303_202'] = Tlvg[2]['Tex']
                outtable[-1]['tau_303_202'] = Tlvg[2]['tau']
                outtable[-1]['Tex_322_221'] = Tlvg[9]['Tex']
                outtable[-1]['tau_322_221'] = Tlvg[9]['tau']
                outtable[-1]['Tex_321_220'] = Tlvg[12]['Tex']
                outtable[-1]['tau_321_220'] = Tlvg[12]['tau']
                outtable[-1]['Temperature'] = T
                outtable[-1]['Column'] = N
                outtable[-1]['nH2'] = n
                outtable[-1]['FWHM'] = dV

outtable.write('ph2cogrid.fits',format='fits',overwrite=True)

