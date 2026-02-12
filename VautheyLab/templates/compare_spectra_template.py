from VautheyLab.standard import *
from VautheyLab.transient_absorption import *
from VautheyLab.miscellaneous import *

# final spectrum
c = Compare_Spectra(files=[''], 
                    labels=[''],
                    colors=[''],
                    scatter=[()],
                    norm=True, 
                    yticks=True,
                    zeroline=True,
                    delays=[0.3], 
                    cuts=[(300, 800)])
c.show()