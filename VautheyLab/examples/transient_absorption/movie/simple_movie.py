from VautheyLab.transient_absorption import Movie

m = Movie(files=['dA1.npy'],
          scatter=[(520, 560)],     # scatter region to exclude (optional)
          experiment='femto',       # experiment type 
          t_cuts=[0.3, 1800],       # time cuts in ps if experiment = 'femto' else in ns (optional)
          wl_cuts=[330, 780],       # wavelength cuts in nm (optinal)
          colors=['crimson'],   
          labels=[r'your molecule'],
          movname='movie.mp4',      # file name to save
          ylim=[-12.5, 25.5],       # file name to save (limits of y-axis, optional)
          before=True,              # decide to plot traces before (optional)
          steady_state=[['abs.txt', (330, 700), -9.5, 'b', 'Abs.']])

# steady state spectra provided as lists [filename, cuts, scale factor, label]
# absorption/emission spectra can be exported from respectrive classes with export=True 

m.render()