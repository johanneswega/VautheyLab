from VautheyLab.transient_absorption import Compare_Overviews

o = Compare_Overviews(files=['dA1.npy', 'dA2.npy'],
                      figsize=(11, 5),      # figure size (optional)
                      t_cuts=[0.3, 1800],   # time cuts
                      titles=['mol. 1', 'mol. 2'], # subplot titles
                      delays = [0.2, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800], # time delays to plot (optional)
                      scatter=[(520, 560), (520, 560)],     # scatter (optional)
                      norm=True)                            # normalize (optional)
o.show()