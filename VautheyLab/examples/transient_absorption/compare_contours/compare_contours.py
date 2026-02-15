from VautheyLab.transient_absorption import Compare_Contours

c = Compare_Contours(files=['dA1.npy', 'dA2.npy'], 
                    figsize=(10, 5),                # figure size (optional)
                    titles=['mol. 1', 'mol. 2'],    # subplot titles
                    t_cuts=[0.3, 1800],   
                    wl_cuts=[340, 750],
                    yscale='log',       # y-axis scale (optional)
                    scale=[-1.2, 1.2],  # scale of the colorbar
                    lines=False,        # whether to plot black lines on contours (optional)
                    nlevels=55,         # number of contours (optional)
                    scatter=[(520, 560), (520, 560)],     # scatter (optional)
                    norm=True,                            # normalize (optional)
                    colorbar=False)                       # whether to plot colorbar (optional)
c.show()