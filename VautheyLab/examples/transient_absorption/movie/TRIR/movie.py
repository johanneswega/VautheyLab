from VautheyLab.transient_absorption import Movie

m = Movie(files=['dA.pdat'],
        t_cuts=[0.3, 500],
        IR=True, # for TRIR
        before=True,  
        figsize=(5, 3.5), # size of the figure (optional)
        labels=[r'TCNQ$^{\bullet –}$'],
        colors=['r'], 
        steady_state=[['FTIR.txt', (1e3, 3e3), -8, 'b', 'FTIR']],
        movname='TCNQ.mp4', 
        ylim=[-10, 4]) # limit of the y-axis plot

m.render() 
