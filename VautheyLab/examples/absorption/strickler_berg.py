from VautheyLab.steady_state import Absorption

# concentration in M
c = 24e-6
# pathlength in cm
l = 1

# initialize absorption class
a = Absorption(files=['abs_file.csv'],
             cuts=[(350, 800)], savefig=True,
             colors=['r'],
             conc=[24e-6], 
             pathlegnths=[1],
             labels=['your molecule'])

# use calc_oscillator_strength(limits, n, nu0, file_index)
# integration limits as list in nm
# refractive index 
# center frequency in kk
a.calc_oscillator_strength([440, 600], 1.421, 18.8, 0)
a.show()