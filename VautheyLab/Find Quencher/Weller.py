import numpy as np

# define type of quenching
print("")
quench = input('    Reductive or Oxidative Quenching? [r / o]: ')

if quench=='o':
    print("")
    E0 = float(input(   'E_00 / eV = '))
    E1 = float(input(   'E0 (C+ / C) = '))
    E2 = float(input(   'E0 (Q / Q-) = '))
    dG = -1*(E2 - E1 + E0)
    print("")
    print("dG0 = %.3g eV"%(dG))

elif quench=='r':
    E0 = float(input(   'E_00 / eV = '))
    E1 = float(input(   'E0 (C / C-) = '))
    E2 = float(input(   'E0 (Q+ / Q) = '))
    dG = -1*(E1 + E0 - E2)
    print("dG0 = %.3g eV"%(dG))
    print("")