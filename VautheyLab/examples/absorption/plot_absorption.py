from VautheyLab.steady_state import Absorption

a = Absorption(files=['abs_file.csv'],
             cuts=[(350, 800)],
             colors=['r'], savefig=True,
             labels=['your molecule'],
             norm=True)
a.show()
