from astropy.table import Table
import pyradex
import numpy as np

R = pyradex.Radex(column=1e16,collider_densities={'pH2':1e4,'oH2':0},species='p-nh3',temperature=15)

fortho = 0.75 # OPR = 1:3

nFWHM = 5
FHWMmin = 0.5
FHWMmx = 5

nDens = 21
nlower = 2
nupper = 6

nCol = 41
Nlower = 10
Nupper = 16

nTemp = 21
Tlower = 10
Tupper = 100

Temps = np.logspace(1,2.5,nTemp)
Cols = 1e1**np.linspace(Nlower,Nupper,nCol)
Densities = 1e1**(np.linspace(nlower,nupper,nDens))
FWHM = np.logspace(np.log10(0.5),np.log10(5),nFWHM)

outtable = Table(names = ['Tex_11','Tex_22','Tex_33','Tex_44',
                          'tau_11','tau_22','tau_33','tau_44',
                          'Temperature','Column','nH2','FWHM'])

for T in Temps:
    for N in Cols:
        for n in Densities:
            for dV in FWHM:
                #opr = 9*np.exp(-170.5/T)
                #fortho = opr/(1+opr)
                fortho = 0
                Tlvg = R(collider_densities={'oH2':0,'pH2':(1-fortho)*n}, column=N,
                         species='p-nh3',temperature=T,deltav=dV)
                outtable.add_row()
                outtable[-1]['Tex_11'] = Tlvg[8]['Tex']
                outtable[-1]['tau_11'] = Tlvg[8]['tau']
                outtable[-1]['Tex_22'] = Tlvg[9]['Tex']
                outtable[-1]['tau_22'] = Tlvg[9]['tau']
                outtable[-1]['Tex_44'] = Tlvg[10]['Tex']
                outtable[-1]['tau_44'] = Tlvg[10]['tau']
                outtable[-1]['Temperature'] = T
                outtable[-1]['Column'] = N
                outtable[-1]['nH2'] = n
                outtable[-1]['FWHM'] = dV
                Tlvg = R(collider_densities={'oH2':0,'pH2':(1-fortho)*n}, column=N,
                         species='o-nh3',temperature=T,deltav=dV)
                outtable[-1]['Tex_33'] = Tlvg[5]['Tex']
                outtable[-1]['tau_33'] = Tlvg[5]['tau']
                
outtable.write('nh3grid.fits',format='fits',overwrite=True)

