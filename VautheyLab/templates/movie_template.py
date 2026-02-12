from VautheyLab.transient_absorption import Movie

m = Movie(files=[''],
          scatter=[()],
          colors=[''],
          labels=[''],
          t_cuts=[0.3, 1800],
          ylim=[-10, 10],
          movname='movie.mp4',
          before=True,
          norm=True,
          yticks=False,
          steady_state=[['', (), 1, 'b', 'Abs.']])
m.render()