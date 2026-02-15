from VautheyLab.transient_absorption import Contour

c = Contour(file='dA.npy', 
            scale=[-27, 27],        # scale of the colorbar
            scatter=[520, 560],     # scatter (optional)
            cmap='RdBu_r',          # color map (optional)
            t_cuts=[-0.3, 1800],
            wl_cuts=[350, 750],
            nlevels=56,             # number of contour levels (optional)
            lines=True,             # whether to draw black lines on the contours (optional)
            ylim=[-0.4, 1800],      # limits of the y-axis scale (optional)
            yscale='symlog')        # set y-scale (optional)
c.show()