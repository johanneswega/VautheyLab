import numpy as np

# give excited state energy and oxidation potential of your chromophore vs. SCE
E00 = 2.85
Eox = 0.85

# load donors
Accpetors = []
Ered_A = []

fh = open('acceptors.txt', 'r')
for line in fh:
    # get redution potential
    if not line[line.find('.')-1:line.find('.')+3]=='':
        # only multiply with -1 with no + in potential
        if not '+' in line:
            Ered_A.append(-1*float(line[line.find('.')-1:line.find('.')+3]))
        else:
            Ered_A.append(float(line[line.find('.')-1:line.find('.')+3]))
        Accpetors.append(line[2:line.find('.')-2])
fh.close()

# calculate driving force in eV
Ered_A = np.array(Ered_A)
dG = (Eox - Ered_A) - E00

# sort lists
sorted_lists = sorted(zip(dG, Accpetors))
dG, Accpetors = zip(*sorted_lists)

# print results
print("")
print("C* + Q --> C+ + Q-")
print("")
print("Chromophore parameters: E00 = %.3g eV / Eox = %.3g vs SCE"%(E00, Eox))
print("")
print("Possible Quenchers/Acceptors for Cation Generation:")
print("")
for i in range(len(dG)):
    if dG[i]<0.3:
        print('dG = %.3g eV     %s'%(dG[i], Accpetors[i]))
print("")