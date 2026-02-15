from VautheyLab.standard import *
from VautheyLab.transient_absorption import *
from VautheyLab.miscellaneous import *

# final spectrum
c = Compare_Spectra(files=['dA1.npy', 'dA2.npy'],
                    labels=['mol. 1', 'mol. 2'],
                    colors=['r', 'b'],
                    delays=[0.3],                       # time delay to plot spectrum
                    scatter=[(520, 560), (520, 560)],   # scatter (optional)
                    norm=True,                          # normalize spectra to max (optional)
                    cuts=[360, 720],                    # wavelength cuts
                    MA=[True, True],                    # use a moving average filter (optional)
                    MA_npoints=[10, 10],                  # how many points to use for MA filter (optional)
                    steady_state=[['abs.txt', (360, 700), -0.5, 'b', 'Abs.']])

# steady state spectra (optional) provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with export=True 
c.show()