from VautheyLab.transient_absorption import Compare_Contours

c = Compare_Contours(files=[],
                     titles=[],
                     t_cuts=[0.3, 100],
                     wl_cuts=[340, 750],
                     yscale='log',
                     scale=[-1, 1],
                     scatter=[(100,200)],
                     colorbar=False)
c.show()