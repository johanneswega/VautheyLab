from VautheyLab.transient_absorption import Global_Analysis

g = Global_Analysis(file='data.pdat',
                    IR=True,                            # set true for TRIR data 
                    model='sequential',                 # choose model
                    IRF=False,                          # whether to do a convolutional fit with analytical gaussian (optional)
                    p0=[1, 10, 100],                    # guess parameters (lifetimes) | if IRF==True p = [lifetimes, t0, FWHM]
                    wavelengths=[2100, 2150, 2200],     # wavenumbers in cm-1 to plot fitted traces (optional)
                    t_cuts=[0.3, 250],                  # limit time axis (optional)
                    wl_cuts=[2100, 2250],               # limit wavelength axis in cm-1 (optional)
                    scale=[-15, 15])                    # scale of contour plots
g.fit()