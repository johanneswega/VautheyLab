from VautheyLab.standard import *

# sellmeier parameters for matrials
# look which sellmeier eqn to use
#
# Quartz: Sellmeier1
# ACN: Sellmeier2
# CHF (chloroform): Sellmeier3 
# CHX: Sellmeier4 
# THF: Sellmeier3
# EtOH: Sellmeier3
# OctOH: Sellmeier3
# Hex: Sellmeier5
# BZN: Sellmeier5
# TOL: Sellmeier5
# BK7: Sellmeier1
# SF10: Sellmeier1
# LAK21: Sellmeier1
# CaF2: Sellmeier1

# types of Sellmeier equations: 
# 1.) n = [A*x**2/(x**2 - B) + C*x**2/(x**2 - D) + E*x**2/(x**2 - F) + 1]**(0.5)
# 2.) n = (A - B*wl**2 + C*wl**(-2) + D*wl**(-4) + E*wl**(-6))**0.5
# 3.) n = ((A*x**2/(x**2 - B)) + (C*x**2/(x**2 - D)) + 1)**(0.5)
# 4.) n = A + B*x**(-2) + C*x**(-4)
# 5.) n = ((A*x**2 / (x**2 - B)) + 1)**(0.5)

material = {'quartz': [0.6961663, 0.0684043**2, 0.4079426, 0.1162414**2, 0.8974794, 9.896161**2],
            'ACN': [1.78207648, 0.001147808, 0.003252916, 0.000655275, 0.000001251],
            'CHF': [1.04647, 0.01048, 0.00345, 0.15207], 
            'CHX': [1.41545, 0.00369, 0.00004],
            'THF': [0.9544, 0.0935**2, 0.0147, 3.450**2],
            'EtOH': [0.8312, 0.0964**2, 0.0192**2, 2.97**2],
            'OctOH': [1.013, 0.0965**2, 0.0112, 2.90**2],
            'Hex': [0.8660, 0.0988**2],
            'BZN': [1.180, 0.138**2],
            'TOL': [1.171, 0.134**2],
            'BK7': [1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653],
            'SF10': [1.62153902, 0.0122241457, 0.256287842, 0.0595736775, 1.64447552, 147.468793],
            'LAK21': [1.22718116, 0.00602075682, 0.420783743, 0.0196862889, 1.01284843, 88.4370099],
            'CaF2': [0.5675888, 0.050263605**2, 0.4710914, 0.1003909**2, 3.8484723, 34.649040**2]}

# sellmeier for Hex
def Sellmeier5(wl, A, B):
    wl = wl/1000
    x = wl 
    n = ((A*x**2 / (x**2 - B)) + 1)**(0.5)
    return n

def Sellmeier5_first_derivative(wl, A, B):
    wl = wl/1000
    x = wl
    dndwl = -(A * B * x) / ((x**2 - B)**2 * np.sqrt((A * x**2) / (x**2 - B) + 1))
    return dndwl

def Sellmeier5_second_derivative(wl, A, B):
    wl = wl/1000
    x = wl
    dn2dwl2 = (A * B * ((3 * A + 3) * x**4 - 2 * B * x**2 - B**2)) / ((x**2 - B)**4 * ((A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

# sellmeier for CHX
def Sellmeier4(wl, A, B, C):
    wl = wl/1000
    x = wl
    n = A + B*x**(-2) + C*x**(-4)
    return n

def Sellmeier4_first_derivative(wl, A, B, C):
    wl = wl/1000
    x = wl
    dndwl = -(2 * B * x**2 + 4 * C) / x**5
    return dndwl

def Sellmeier4_second_derivative(wl, A, B, C):
    wl = wl/1000
    x = wl
    dn2dwl2 = (6 * B * x**2 + 20 * C) / x**6
    return dn2dwl2

# type of Sellmeier of chloroform 
def Sellmeier3(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    n = ((A*x**2/(x**2 - B)) + (C*x**2/(x**2 - D)) + 1)**(0.5)
    return n

def Sellmeier3_first_derivative(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    dndwl = ((2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2) / (2 * np.sqrt((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1))
    return dndwl

def Sellmeier3_second_derivative(wl, A, B, C, D):
    wl = wl/1000
    x = wl
    dn2dwl2 = ((2 * C) / (x**2 - D) - (10 * C * x**2) / (x**2 - D)**2 + (8 * C * x**4) / (x**2 - D)**3 + (2 * A) / (x**2 - B) - (10 * A * x**2) / (x**2 - B)**2 + (8 * A * x**4) / (x**2 - B)**3) / (2 * np.sqrt((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)) - ((2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2)**2 / (4 * ((C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

# type of Sellmeier used for ACN
def Sellmeier2(wl, A, B, C, D, E):
    wl = wl/1000
    n = (A - B*wl**2 + C*wl**(-2) + D*wl**(-4) + E*wl**(-6))**0.5
    return n

def Sellmeier2_first_derivative(wl, A, B, C, D, E):
    wl = wl / 1000  # Convert to microns
    x = wl
    dn_dwl = -(B * x**8 + C * x**4 + 2 * D * x**2 + 3 * E) / (x**7 * np.sqrt(-B * x**2 + C / x**2 + D / x**4 + E / x**6 + A))
    return dn_dwl

def Sellmeier2_second_derivative(wl, A, B, C, D, E):
    wl = wl / 1000  # Convert to microns
    x = wl
    d2n_dwl2 = -(A * B * x**14 + 6 * B * C * x**12 + (15 * B * D - 3 * A * C) * x**10 + (28 * B * E - 10 * A * D - 2 * C**2) * x**8 + (-21 * A * E - 9 * C * D) * x**6 + (-18 * C * E - 6 * D**2) * x**4 - 19 * D * E * x**2 - 12 * E**2) / (x**14 * (-B * x**2 + C / x**2 + D / x**4 + E / x**6 + A)**(3 / 2))
    return d2n_dwl2

# type of Sellmeier equation for quartz
def Sellmeier1(wl, A, B, C, D, E, F):
    # wavelength is in um
    wl = wl/1000
    x = wl
    n = (A*x**2/(x**2 - B) + C*x**2/(x**2 - D) + E*x**2/(x**2 - F) + 1)**(0.5)
    return n

# in um^-1
def Sellmeier1_first_derivative(wl, A, B, C, D, E, F):
    wl = wl/1000
    x = wl
    dndwl = ((2 * E * x) / (x**2 - F) - (2 * E * x**3) / (x**2 - F)**2 + (2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2) / (2 * np.sqrt((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1))
    return dndwl

# in um^-2
def Sellmeier1_second_derivative(wl, A, B, C, D, E, F):
    wl = wl/1000
    x = wl
    dn2dwl2 = ((2 * E) / (x**2 - F) - (10 * E * x**2) / (x**2 - F)**2 + (8 * E * x**4) / (x**2 - F)**3 + (2 * C) / (x**2 - D) - (10 * C * x**2) / (x**2 - D)**2 + (8 * C * x**4) / (x**2 - D)**3 + (2 * A) / (x**2 - B) - (10 * A * x**2) / (x**2 - B)**2 + (8 * A * x**4) / (x**2 - B)**3) / (2 * np.sqrt((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)) - ((2 * E * x) / (x**2 - F) - (2 * E * x**3) / (x**2 - F)**2 + (2 * C * x) / (x**2 - D) - (2 * C * x**3) / (x**2 - D)**2 + (2 * A * x) / (x**2 - B) - (2 * A * x**3) / (x**2 - B)**2)**2 / (4 * ((E * x**2) / (x**2 - F) + (C * x**2) / (x**2 - D) + (A * x**2) / (x**2 - B) + 1)**(3 / 2))
    return dn2dwl2

def Group_velocity(wl, n, dndwl):
    wl = wl/1000
    A = sc.c/n
    B = 1 - (wl/n)*dndwl
    vg = A/B
    return vg

def Group_velocity_dispersion(wl, dn2dwl2):
    wl = wl/1000
    # speed of light in um/fs
    cfs = sc.c*(10**6 / (10**(15)))
    # GVD in fs^2 / um
    gvd = (wl**3/(2*np.pi*cfs**2))*dn2dwl2
    # in fs^2 / mm
    gvd = gvd*1000
    return gvd
