from VautheyLab.transient_absorption import Global_Analysis

# suppose we want to model as branching scheme like this one
#                       A
#                [k0]|     | [k1]
#                    v     v
#                    B     C
#                [k2]|     | [k3]
#                    v     v

# the rate equations for the species are 
# dA/dt = -(k0 + k1) * A    + 0 * B         + 0 * C
# dB/dt = +k0 * A           -k2 * B         + 0 * C
# dB/dt = +k1 * A           + 0 * B         - k3 * C

# which can be written in matrix form
#       [A]    [-(k0 + k1),     0,      0] [A]
# d/dt  [B] =  [+k0,          -k2,      0] [B]
#       [C]    [+k1,            0,    -k3] [C]

# for the target analysis we need to provide this K-matrix 
# as a numpy array to the class
import numpy as np
K = np.array([['-k0 -k1', '0', '0'],
              ['+k0', '-k2', '0'],
              ['+k1', '0', '-k3']])

# as well as the starting concentrations at t=0
# assuming everything is in A at t=0
C0 = np.array([1, 0, 0])

g = Global_Analysis(file='dA.npy',
                    model='target',              # choose model
                    K=K,                         # provide rate contant matrix
                    C0=C0,                       # provide inital conditions
                    IRF=False,                   # whether to do a convolutional fit with analytical gaussian (optional)
                    p0=[10, 10, 100, 100],       # guess parameters (lifetimes) | if IRF==True p = [lifetimes, t0, FWHM]
                    wavelengths=[400, 450, 500,
                                 550, 600, 650], # wavelength to plot fitted traces (optional)
                    t_cuts=[5, 1800],            # limit time axis (optional)
                    wl_cuts=[350, 750],          # limit wavelength axis (optional)
                    scatter=[520, 560],          # exclude scatter (optional)
                    scale=[-25, 25])             # scale of contour plots

g.fit()