import numpy as np

# input 
print("")
inp = input('Laser Quantatiy: ')
# check units 
d = float(input('Beam diameter [um]: '))*10**-4 
r = d/2
A = np.pi*r**2

if 'J' in inp and not 'cm2' in inp: 
    r = float(input('Rep-rate [kHz]: '))
    # calculate in Hz 
    r = r*1000
    # calculate energy per pulse
    if 'u' in inp:
        E = float(inp[:inp.find('uJ')])
        E = E*10**-6
        P = E*r
    elif 'm' in inp:
        E = float(inp[:inp.find('mJ')])
        E = E*10**-3
        P = E*r   
    else:
        E = float(inp[:inp.find('J')])
        P = E*r        
    print("")   
    if P>=10**(-3) and P<1:
        print("E = %.3g J per pulse --> P = %.3g mW at %.3g kHz"%(E, P*10**3, r/1000))
    elif P>=10**(-6) and P<10**(-3):
        print("E = %.3g J per pulse --> P = %.3g uW at %.3g kHz"%(E, P*10**6, r/1000))
    else:
        print("E = %.3g J per pulse --> P = %.3g W at %.3g kHz"%(E, P, r/1000))
    print("")
    # calculate fluence in uJ/cm^2
    F = E*10**6/A
    print("Fluence F = %.3g uJ/cm^2 = %.3g mJ/cm^2"%(F, F*10**-3))

if 'W' in inp: 
    r = float(input('Rep-rate [kHz]: '))
    # calculate in Hz 
    r = r*1000
    # calculate energy per pulse
    if 'u' in inp:
        P = float(inp[:inp.find('uW')])
        P = P*10**-6
        E = P/r
    elif 'm' in inp:
        P = float(inp[:inp.find('mW')])
        P = P*10**-3
        E = P/r        
    else:
        P = float(inp[:inp.find('W')])
        E = P/r        
    print("")   
    if E>=10**(-3) and E<1:
        print("P = %.3g W at %.3g kHz --> E = %.3g mJ per pulse"%(P, r/1000, E*10**3))
    elif E>=10**(-6) and E<10**(-3):
        print("P = %.3g W at %.3g kHz --> E = %.3g µJ per pulse"%(P, r/1000, E*10**6))
    else:
        print("P = %.3g W at %.3g kHz --> E = %.3g J per pulse"%(P, r/1000, E))
    print("")
    # calculate fluence in uJ/cm^2
    F = E*10**6/A
    print("Fluence F = %.3g uJ/cm^2 = %.3g mJ/cm^2"%(F, F*10**-3))

if 'cm2' in inp:
    r = float(input('Rep-rate [kHz]: '))
    # calculate in Hz 
    r = r*1000
    # calculate E 
    if 'uJ' in inp:
        F = float(inp[:inp.find('u')])
        E = A*F/10**6
    if 'mJ' in inp:
        F = float(inp[:inp.find('m')])
        E = A*F/10**3
    P = E*r
    print('''E = %.3g J = %.3g mJ = %.3g uJ per pulse 
          --> P = %.3g W = %.3g mW = %.3g uW at %.3g kHz'''%(E, E*10**3, E*10**6, P, P*10**3, P*10**6, r))
    
print("")