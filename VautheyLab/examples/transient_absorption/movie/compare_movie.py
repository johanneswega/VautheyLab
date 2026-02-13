from VautheyLab.transient_absorption import Movie

m = Movie(files=['NP.npy', 'Ac.npy'],
          scatter=[(520, 560), (520, 560)],
          colors=['darkorange', 'crimson'],
          labels=[r'R = $n$-Pr', r'R = Az'],
          movname='comp.mp4',
          t_cuts=[0.3, 1800],
          wl_cuts=[330, 750],
          ylim=[-0.65, 1.2],
          before=False,
          norm=True, # normalize spectra (optional)
          normat=[0.2, 380], # noramlize at specific (delay, wavelength)
          yticks=False) # remove y-ticks of the plot

m.render()