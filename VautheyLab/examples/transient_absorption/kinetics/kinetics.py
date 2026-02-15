from VautheyLab.transient_absorption import Kinetics

k = Kinetics(file='dA.npy',
            wavelengths=[400, 500, 600], # wavelength to plot kinetics
            outside=True, # put legend next to plot (optional)
            colors=['b', 'g', 'r'], 
            cuts=[0.3, 1800], # time cuts
            xscale='log') # set x-axis scale to logarithmic
k.show()