from VautheyLab.standard import *

# function to convert nm to kk and vice versa
def nm_to_kk(val):
    return (1/val)*10**4

# function to convert cm-1 to fs^-1
def cm_to_fs(val):
    return val*2.998e+10*10**(-15)

# function to convert cm-1 to ps^-1
def cm_to_ps(val):
    return (sc.c * (val*100))*10**(-12)

# FT from ps to cm-1
def FT_ps_to_cm(t, E):
    # Do the FFT
    fhat = np.fft.fftshift(np.fft.fft(E))
    freq = np.fft.fftshift(np.fft.fftfreq(len(t), d=t[1]-t[0]))
    f = freq[freq > 0]
    wn = (1 / (sc.c / f * 10**-12)) * 10**-2
    FT = np.abs(fhat[freq > 0])
    return wn, FT

# print fit parameters with erros
def print_params(p, pcov):
    print('')
    print('### Optimized Fit Parameters ###')
    print('')
    for i in range(len(p)):
        err = pcov[i,i]**(0.5)
        print('%i. %.3g +- %.3g'%(i+1, p[i], err))
    print('')

def sort_with_respect(sort_after, *args):
    # Sort all lists based on sort_after
    sorted_lists = sorted(zip(sort_after, *args))  
    
    # Unpack sorted tuples
    sorted_arrays = list(zip(*sorted_lists))  

    # Convert to numpy arrays and return
    return tuple(np.array(arr) for arr in sorted_arrays)
    
# Moving Average Filter
def moving_average(x, window_size):
    window = np.ones(window_size) / window_size
    return np.convolve(x, window, mode='valid')

# Function to calculate rainbow list
def rainbow(liste, r=None):
    if r!=None:
        return plt.cm.rainbow_r(np.linspace(0,1,len(liste)))
    else:
        return plt.cm.rainbow(np.linspace(0,1,len(liste)))
    
def gist_rainbow(liste, r=None):
    if r!=None:
        return plt.cm.gist_rainbow_r(np.linspace(0,1,len(liste)))
    else:
        return plt.cm.gist_rainbow(np.linspace(0,1,len(liste)))
    
# function to find index
def find_index(x, find):
    return np.argmin(np.abs(x - find))

# function for reading theoretical IR spectrum
def get_peaks(fname):
    x = []
    y = []
    with open(fname, "r") as file: 
        for line in file:
            # only read lines starting with #
            if line.startswith('#'):
                # split line
                line = line.split()
                # if second entry can be converted to a number add it
                try:
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                except ValueError:
                    continue
    return np.array(x), np.array(y)

# load theoretical IR
def load_theo_IR(fname, cuts, fwhm, scale):
    # load peak information
    wnp, Ip = get_peaks(fname)

    wnp = scale*wnp

    # take only those peaks that lie within the range of cuts of exp. 
    Ip = Ip[(wnp>cuts[0]) & (wnp<cuts[1])]
    wnp = wnp[(wnp>cuts[0]) & (wnp<cuts[1])]
    
    # make frequency axis
    wn = np.linspace(cuts[0], cuts[1], 5000)
    I = np.zeros(len(wn))

    # construct spectrum
    from VautheyLab.fit_functions import lorentzian

    for i in range(len(Ip)):
        I += lorentzian(wn, Ip[i], wnp[i], fwhm)
    return (1/wn)*10**4, wn, I

# function to export theoIR to .txt as comp steady state in scripts
def export_theoIR(fname, cuts, fwhm, scale):
    head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
    wl, wn, A = load_theo_IR(fname, cuts, fwhm, scale)  
    np.savetxt('TRIR_theoretical.txt', 
                np.column_stack([wl, wn, A]),
                header=head, delimiter=',')  

# calculate theoretical FTIR spectrum 
def calc_theo_TRIR(fGS, fES, cuts, fwhm, scale):
    head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
    wl, wn, AGS = load_theo_IR(fGS, cuts, fwhm, scale)
    wl, wn, AES = load_theo_IR(fES, cuts, fwhm, scale)
    dA = AES - AGS
    np.savetxt('TRIR_theoretical.txt', 
                np.column_stack([wl, wn, dA]),
                header=head, delimiter=',')
    
# get sorted file after integer i.e. 2 before 10
def get_sorted_files(path):
    # Get all file names
    files = os.listdir(path)

    # get rid of .DS_store if present
    files = [i for i in files if not '.DS' in i]

    # Sort numerically (natural order)
    files_sorted = sorted(files, key=lambda x: int(os.path.splitext(x)[0]))
    return np.array([path + '/' + i for i in files_sorted])

# function to convert shida csv to txt to be used in absorption script
def shida_to_txt(shida, ma=False, units='wn'):
    data = np.loadtxt(shida, delimiter=',')
    if units=='wn':
        wn = data[:, 0]
        wl = nm_to_kk(wn)
    else:
        wl = data[:, 0]
        wn = nm_to_kk(wl)        
    A = data[:, 1]
    if ma!=False:
        wl = moving_average(wl, ma)
        wn = moving_average(wn, ma)
        A = moving_average(A, ma)
    # sort with respect to wavelengths
    #A = np.flip(A)
    #wn = np.flip(wn)
    wl, wn, A = sort_with_respect(wl, wn, A)
    np.savetxt('%s.txt'%(shida[:shida.find('.')]), 
                np.column_stack([wl, wn, A]),
                header='wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance', delimiter=',')
    
def get_data(file, delimiter=',', skiprows=1, usecols=[0,1], xcuts=None):
    data = np.loadtxt(file, delimiter=delimiter, skiprows=skiprows, usecols=usecols)
    x = data[:,0]
    y = data[:,1]
    x, y = sort_with_respect(x, y)
    if xcuts!=None:
        y = y[(x>xcuts[0]) & (x<xcuts[1])]
        x = x[(x>xcuts[0]) & (x<xcuts[1])]
    return x, y

def yline(ax, y, color='k', linestyle='--', alpha=0.3):
    ax.axhline(y=y, color=color, linestyle=linestyle, alpha=alpha)

def xline(ax, x, color='k', linestyle='--', alpha=0.3):
    ax.axvline(x=x, color=color, linestyle=linestyle, alpha=alpha)

# solvent dictonary 
# 0 - n
# 1 - eps
# 2 - df
# 3 - eta
# 4 - alpha = Kamlet-Taft hydrogen donor abilitiy

def onsager(n, e_r):
    return (2*(e_r - 1)/(2*e_r+1)) - (2*(n**2 - 1)/(2*n**2+1))

#print(onsager(1.53269 , 12.8))

solvs = {'ACN': [1.341, 35.94, 0.6115, 0.341, 0.19],
         'CHX': [1.4235, 2.02, -0.001497, 0.898, np.nan],
         'BCN': [1.525, 25.2, 0.4725, 1.237, np.nan],
         'TOL': [1.4941, 2.38, 0.02815, 0.553, 0.0],
         'DCM': [1.421, 8.93, 0.4364, 0.411, 0.13],
         'THF': [1.405, 7.58, 0.4207, 0.462, 0.0],
         'DMSO': [1.477, 46.45, 0.5274, 1.991, 0.0],
         'CHF': [1.442, 4.89, 0.3032, 0.536, 0.20],
         'MeOH': [1.3265, 32.66, 0.6186, 0.551, 0.98],
         'EtOH': [1.3594, 24.55, 0.579, 1.083, np.nan],
         'GLY': [1.473, 42.5, 0.527, 945, 1.21],
         'iBuOH': [1.3939, 17.93, 0.5326, 3.333, 0.79],
         'BuOH': [1.3974, 17.51, 0.5283, 2.573, 0.84],
         'DecOH': [1.435, 8.1, 0.4417, 11.32, 0.70],
         'NonOH': [1.4338 , 9.0, 0.4290, 10.27, 0.75],
         'ANIS': [1.5143, 4.33, 0.2265, 0.984, np.nan],
         'PhEtOH': [1.53269 , 12.8, 0.4137, 7.58, np.nan],
         'FUR': [1.4187, 2.94, 0.1609, 0.361, np.nan],
         'MNP': [1.6120, 2.70, onsager(1.6120, 2.70), 2.90, np.nan],
         'BNZ': [1.4979, 2.27, 0.00519, 0.603, np.nan],
         'pXYL': [1.4933, 2.27, 0.007945, 0.605, np.nan],
         'BRF': [1.595, 4.39, 0.186, 1.868, np.nan],
         'H2O': [1.3325, 78.36, 0.6402, 0.8903, 1.17]}