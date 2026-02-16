from VautheyLab.transient_absorption import Global_Analysis

g = Global_Analysis(file='dA.npy',
                    model='sequential',          # choose model
                    IRF=False,                   # whether to do a convolutional fit with analytical gaussian (optional)
                    p0=[1, 10, 100],             # guess parameters (lifetimes) | if IRF==True p = [lifetimes, t0, FWHM]
                    wavelengths=[400, 450, 500,
                                 550, 600, 650], # wavelength to plot fitted traces (optional)
                    t_cuts=[5, 1800],            # limit time axis (optional)
                    wl_cuts=[350, 750],          # limit wavelength axis (optional)
                    scatter=[520, 560],          # exclude scatter (optional)
                    scale=[-25, 25])             # scale of contour plots

g.fit()