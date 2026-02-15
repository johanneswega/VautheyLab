from VautheyLab.transient_absorption import Overview

o = Overview(file='dA.npy',
             delays = [0.2, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800], # time delays to plot
             scatter=[520, 560],       # scatter region to exclude (optional)
             experiment='femto',       # experiment type (optional)
             IR=False,                 # not TRIR data (optional)
             cuts=[330, 780],          # wavelength cuts (optional)
             steady_state=[['abs.txt', (330, 700), -9.5, 'b', 'Abs.']])

# steady state spectra (optional) provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with export=True 

# if the delays argument is not given 
# it is set to: 
# delays = [-1, 0.2, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800] 
# for femto
# and:
# [2, 5, 10, 20, 50, 80, 100, 250, 200, 500, 1e+3, 5e+3, 20e+3, 50e+3, 100e+3, 200e+3, 500e+3] 
# for nano
# by standard

o.show()