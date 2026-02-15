from VautheyLab.transient_absorption import Compare_Kinetics

k = Compare_Kinetics(files=['dA1.npy', 'dA2.npy'],
                     labels=['mol. 1', 'mol. 2'],
                     colors=['b', 'r'],
                     wavelengths=[400, 500],          # wavelength to plot kinetics at
                     cuts=[(0.3, 1800), (0.3, 1800)], # wavelength cuts (optional)
                     xscale='log')                    # scale of the x-axis
k.show()