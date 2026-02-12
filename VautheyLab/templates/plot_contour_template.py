from VautheyLab.transient_absorption import Contour

c = Contour(file='',
            scale=[-5,5],
            scatter=[],
            t_cuts=[-0.2, 1800],
            wl_cuts=[350, 750],
            yscale='symlog')
c.show()