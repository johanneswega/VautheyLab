import numpy as np
from VautheyLab.miscellaneous import *

### ultimate concentration calculator

### c in mM
### rho in g/mL
### V in mL
### M in g/mol
### eps in M^-1 cm^-1
### l in cm
### m in g
### tau in ns

# function to print units 
def print_units():
    units = '''### c in mM
### rho in g/mL
### V in mL
### M in g/mol
### eps in M^-1 cm^-1
### l in cm
### m in g
### tau in ns'''
    print(units)

# function to calculate dilution of a stock solution to 
# a new solution with different concentation
def dilute(cS, cwant, Vcuv, show=True):
    print("")
    Vadd = -cwant*Vcuv/(cwant - cS)
    Vtotal = Vadd + Vcuv
    if show==True:
        print("""To make a c_want = %.3g mM solution from c_stock = %.3g mM stock solution. Take V_cuv = %.3g mL
            of pure solvent and add Vadd = %.3g mL = %.3g uL of stock solution. The total solution volume will be %.3g mL"""
            %(cwant, cS, Vcuv, Vadd, Vadd*1000, Vtotal))
        print("")
    return Vadd*1000

# function to calculate new concentration after dilution
def conc_after_dil(cS, Vadd, Vcuv):
    print("")
    cdil = (cS*Vcuv)/(Vadd + Vcuv)
    Vtotal = Vadd + Vcuv  
    print("""After adding Vadd = %.3g mL of solvent to a stock solution with cS = %.3g mM and Vcuv = %.3g mL the new solution 
          has a concentration of cdil = %.3g mM and a total volume of Vtot = %.3g mL."""
        %(Vadd, cS, Vcuv, cdil, Vtotal))
    print("")  
    return cdil

# function to calculate new concentration after dilution
def conc_after_dil2(cS, Vadd, Vcuv):
    print("")
    cdil = (cS*Vadd)/(Vadd + Vcuv)
    Vtotal = Vadd + Vcuv  
    print("""After adding Vadd = %.3g mL of stock solution with cS = %.3g mM to a cuvette with Vcuv = %.3g mL the new solution 
          has a concentration of cdil = %.3g mM = %.3g µM and a total volume of Vtot = %.3g mL."""
        %(Vadd, cS, Vcuv, cdil, cdil*1000, Vtotal))
    print("")  
    return cdil

# function to calculate needed concentration to have 
# a certain absobance
def conc_for_abs(A, eps, l):
    print("")
    c = (A/(eps * l))*1000
    print("""For A = %.3g in a l = %.3g cm cuvette you need a c = %.3g mM = %.3g uM 
          if eps = %.3g M-1 cm-1 at this wavelength."""
          %(A, l, c, c*1000, eps))
    print("")
    return c

# calculate concentration from absorbance
def abs_for_conc(c, eps, l):
    print("")
    A = c*10**-3*eps*l
    print("For c = %.3g mM, you need A = %.3g in l = %.3g cm if eps = %.3g M-1 cm-1."%(c, A, l, eps))
    print("")

# function to calculate needed mass of solid substance 
# to make certain concentration
def mass_for_conc(c, V, M):
    m = c*10**(-3) * V* 10**(-3) * M
    print("")
    print("""To make a c = %.3g mM solution, dissolve m = %.3g g = %.3g mg of substance
          in V = %.3g mL of solvent."""%(c, m, m*1000, V))
    print("")

# function to calculate concentration from mass
def conc_from_mass(m, V, M):
    c = (m/(M*V*10**-3))*1000
    print("")
    print("""m = %.3g mg of substance dissolved in V = %.3g mL of solvent corresponds 
            to a concentration of c = %.3g mM"""%(m*1000, V, c))
    print("")
    return c

# function to calculate required volume for certrain conc. and mass
def vol_from_mass(m, c, M):
    V = (m/(M*c*10**(-3)))*1000
    print("")
    print("""For m = %.3g mg of substance dissolved in V = %.3g mL of solvent corresponds 
            to a concentration of c = %.3g mM"""%(m*1000, V, c))
    print("")

# function to calculate how much liquid quencher to add
def liq_quench(cQ, Vcuv, MQ, rhoQ):
    # convert to g/L
    rhoQ = rhoQ * 1000
    # convert to mol/L
    cQ = cQ/1000
    # convert to L
    Vcuv = Vcuv/1000
    # calculate Vadd in L 
    Vadd = -1*cQ*MQ*Vcuv/(cQ*10**(-3)*MQ - rhoQ)
    # convert to mL 
    Vadd = Vadd*1000
    print("")
    print("""If you have Vcuv = %.3g mL in your cuvette 
          and you want cQ = %.3g mM. Directly pipette Vadd = %.3g uL
          of liquid quencher into your cuvette."""%(Vcuv*1000, cQ*1000, Vadd*1000))
    print("")
    return Vadd

# calculate stock solution of pure liquid
def pure_liq(rhoQ, MQ):
    # calculate "stock" concentration of liquid quencher
    cSQ = (rhoQ / MQ)*10**6
    print("Pure liq. quencher has a molar (stock) concentration of %.3g M"%(cSQ/1000))
    print("")
    return cSQ

# function for a titration experiment
def titration(cQ, Vcuv, MQ=None, cSQ=None, rhoQ=None):
    print("")
    # calculate amounts to add each time
    cadd = [cQ[0]]
    for i in range(1, len(cQ)):
        cadd.append(cQ[i] - np.sum(cadd))
    # for liq. quencher
    if cSQ==None:
        # calculate "stock" concentration of liquid quencher
        cSQ = (rhoQ / MQ)*10**6
        print("Pure liq. quencher has a molar (stock) concentration of %.3g M"%(cSQ/1000))
        print("")
    # loop through concentrations
    Vnew = Vcuv
    for i in range(len(cadd)):
        V_add = dilute(cSQ, cadd[i], Vcuv, show=False) 
        Vnew += V_add/1000
        print("%i. For cQ = %.3g mM add V_add = %.3g uL    (V_total = %.3g mL)"%(i+1, cQ[i], V_add, Vnew))
    print("")

# function to calculate concentration for equal small amounts 
# of quencher added
def titra_conc_from_vol_eql(csQ, Vadd, Vcuv, n):
    ctot = 0
    cQ = []
    print("")
    print("""Adding Vadd = %.3g uL of Quencher with cSQ = %.3g mM to a solution with Vcuv = %.3g mL each time
           results in the following concentrations:"""%(Vadd*1000, csQ, Vcuv))
    print("")
    for i in range(n):
        Vcuv += Vadd
        c_each = csQ*Vadd/(Vcuv)
        ctot += c_each
        print("%i. Addition [Q] = %.3g mM = %.3g uM"%(i+1, ctot, ctot*1000))
        cQ.append(ctot*1000)
    print("")
    print(cQ)
    return cQ

# function calculate concentration for arbitrary amounts of quencher stock added to cuvette
def titra_conc_from_vol(csQ, Vadd, Vcuv):
    ctot = 0
    cQ = []
    print("")
    print("""Adding Vadd = x uL of Quencher with cSQ = %.3g mM to a solution with Vcuv = %.3g mL
           results in the following concentrations:"""%(csQ, Vcuv))
    print("")
    for i in range(len(Vadd)):
        Vcuv += Vadd[i]
        cnow = csQ*Vadd[i]/(Vcuv) 
        ctot += cnow
        print("%i. Addition Vadd = %.3g uL => [Q] = %.3g mM = %.3g uM        | new solution volume = %.3g mL"%(i+1, Vadd[i]*1000, ctot, ctot*1000, Vcuv))
        cQ.append(ctot*1000)
    print("")
    return cQ        

# function to calculate what quencher concentration according to Stern-Volmer
# and diffusion controlled quenching is needed to reduce the lifetime from tau0 to tau
def SV_quench_conc(tau0, tau, solv, T=298.15):
    # get viscosity of solvent
    eta = solvs[solv][3]
    # calculate kdiff
    kd = (8*sc.R*1000*T)/(3*eta/1000)
    # calculate [Q]
    Q = ((tau0/tau) - 1) / (kd * tau0 * 10**-9)
    print('')
    print('k_diff = %.3g M-1 s-1 in %s at T = %.5g K'%(kd, solv, T))
    print('''To reduce the lifetime from tau0 = %.3g ns to tau = %.3g ns you need a quencher 
          concentration of [Q] = %.3g mM = %.3g uM'''%(tau0, tau, Q*10**3, Q*10**6))
    print('')
    return Q*10**6

# function to calculate expected lifetime reduction for diffusion controlled quenching
def SV_tau_from_c(tau0, Q, solv, T=298.15):
    # get viscosity of solvent
    eta = solvs[solv][3]
    # calculate kdiff
    kd = (8*sc.R*1000*T)/(3*eta/1000)
    # calculate tau
    # convert to M
    Q = Q/1000
    tau = tau0 / (1 + kd * tau0 * 10**-9 * Q)
    print('')
    print('k_diff = %.3g M-1 s-1 in %s at T = %.5g K'%(kd, solv, T))
    print('''For a quencher conc. of [Q] = %.3g mM and a lifetime of tau0 = %.3g ns we expect the lifetime
          to reduce to tau = %.3g ns for pure diffusional quenching.'''%(Q*1000, tau0, tau))
    print('')

# function to calculate needed amount of solvent to add to a concentrated 
# stock solution to make wanted concentration
def dilution_series(cstock, Vcuv, cwant):
    print("")
    added = []
    for i in range(len(cwant)):
        # volume to add
        Va = (cstock*Vcuv - cwant[i]*Vcuv)/cwant[i]
        if i==0:
            print("For c = %.3g mM add V = %.3g mL to stock solution of %.3g mM in %.3g mL."%(cwant[i], Va, cstock, Vcuv))
            print("")
        else:
            print("To further dilute to c = %.3g mM add again %.3g mL of solvent. The total volume is now %.3g mL."%(cwant[i], Va - np.sum(added), Vcuv + Va))
            print("")
        added.append(Va)
    print("")

### DMA ###
# M = 121.19 g/mol
# rho = 0.956 g/mL
    
### DCB ###
# M = 128.134 g/mol
    
### MA ###
# M = 98.06 g/mol