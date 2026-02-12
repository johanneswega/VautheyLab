from VautheyLab.standard import *
from VautheyLab.miscellaneous import *

# load donors
Donors = []
Eox_D = []

fh = open('donors.txt', 'r')
for line in fh:
    # get redution potential
    Eox_D.append(float(line[line.find('+')+1:line.find('+')+5]))
    Donors.append(line[2:line.find('+')])
fh.close()

fig, ax = plt.subplots(1,1, figsize=(15, 8))

Eox_D, Donors = sort_with_respect(Eox_D, Donors)
colors = rainbow(Donors)

for i in range(len(Donors)):
    ax.plot(i, Eox_D[i], 'o', alpha=0.6, color=colors[i])
    ax.text(x=i, y=Eox_D[i]+0.1, s=Donors[i], size=8, rotation=90)

#ticks = [i for i in range(len(Donors))]
ax.set_xticks([])
ax.set_ylabel(r'$E^0_{\text{ox}}$ / V vs. SCE')
ax.set_xlabel(r'$\leftarrow$ Electron Donating Ability')
plt.savefig('Donors.svg')
plt.show()
