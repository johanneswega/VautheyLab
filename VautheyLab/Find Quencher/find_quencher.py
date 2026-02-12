import numpy as np

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

# sort with respect to strength
from VautheyLab.miscellaneous import sort_with_respect

Ered_A, Accpetors = sort_with_respect(Ered_A, Accpetors)
Eox_D, Donors = sort_with_respect(Eox_D, Donors)

#print(Ered_A)

print('')
q1 = input('     donor [d] or acceptor [a] ?  ')
print('')

if q1 == 'a':
    Er = input('     E0(A/A-) vs. SCE range = start, stop?   ')
    Ersplit = Er.split(',')
    lim1 = float(Ersplit[0])
    lim2 = float(Ersplit[1])
    
    req = Accpetors[(Ered_A > lim1) & (Ered_A < lim2)]
    Ereq = Ered_A[(Ered_A > lim1) & (Ered_A < lim2)]
    
    print('\n     Molecule Name                          :   E0(A/A-) vs. SCE / V')
    print('      --------------------------------------:----------')
    
    for i in range(len(req)):
        print(f'     {req[i]:<60} : {Ereq[i]:>6} V')

elif q1 == 'd':
    Er = input('     E0(D+/D) vs. SCE range = start, stop?   ')
    Ersplit = Er.split(',')
    lim1 = float(Ersplit[0])
    lim2 = float(Ersplit[1])
    
    req = Donors[(Eox_D > lim1) & (Eox_D < lim2)]
    Ereq = Eox_D[(Eox_D > lim1) & (Eox_D < lim2)]
    
    print('\n     Molecule Name                          :   E0(D+/D) vs. SCE / V')
    print('      --------------------------------------:----------')
    
    for i in range(len(req)):
        print(f'     {req[i]:<60} : {Ereq[i]:>6} V')
else:
    print('input not valid. Only type a or d.')
    quit()