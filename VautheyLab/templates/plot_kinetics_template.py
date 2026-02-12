from VautheyLab.transient_absorption import Kinetics

k = Kinetics(file='',
            ylim=[-2, 2],
            wavelengths=[500],
            colors=['r'],
            cuts=[0.3, 1800],
            xscale='log')
k.show()