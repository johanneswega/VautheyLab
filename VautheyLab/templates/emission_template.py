from VautheyLab.steady_state import Emission

e = Emission(files_em=[''],
             cuts_em=[()],
             labels_em=[''],
             colors_em=['r'],
             files_abs=[''],
             cuts_abs=[()],
             labels_abs=[''],
             colors_abs=['b'],
             norm_abs=True,
             norm_em=True,
             yticks=False,
             corr=True)

e.show()