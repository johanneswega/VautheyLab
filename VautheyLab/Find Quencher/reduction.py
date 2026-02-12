import numpy as np

# give excited state energy and oxidation potential of your chromophore vs. SCE
E00 = 2.4
Ered = -1.29 + 0.38

# load donors
Donors = []
Eox_D = []

fh = open('donors.txt', 'r')
for line in fh:
    # get redution potential
    Eox_D.append(float(line[line.find('+')+1:line.find('+')+5]))
    Donors.append(line[2:line.find('+')])
fh.close()

# calculate driving force in eV
Eox_D = np.array(Eox_D)
dG = (Eox_D - Ered) - E00

# sort lists
sorted_lists = sorted(zip(dG, Donors))
dG, Donors = zip(*sorted_lists)

# print results
print("")
print("C* + Q --> C- + Q+")
print("")
print("Chromophore parameters: E00 = %.3g eV / Ered = %.3g vs SCE"%(E00, Ered))
print("")
print("Possible Quenchers/Donors for Anion Generation:")
print("")
for i in range(len(dG)):
    if dG[i]<1.0:
        print('dG = %.3g eV     %s'%(dG[i], Donors[i]))
print("")