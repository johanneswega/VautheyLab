import numpy as np
from scipy.special import erf
from scipy.special import wofz

# chirp correction function
def koboyashi(l, a, b, c):
    return a + 10**5 * b * l**-2 + 10**6 * c * l**-4

# displaced harmonic oscillator for absorption bands
def DHO(wn, I, S, nu0, om, sigma):
    A = np.zeros(len(wn))
    for m in range(50):
        add = wn*(S**m * np.exp(-S))/(np.math.factorial(m))*np.exp((-(nu0 + (m*om) - wn)**2)/(2*sigma**2))
        A = np.add(A, add)
    return np.array(I*A/np.max(A), dtype=float)
    
# sum of two DHOs
def two_DHOs(wn, A1, S1, nu01, om1, sigma1, A2, S2, nu02, om2, sigma2):
    return DHO(wn, A1, S1, nu01, om1, sigma1) + DHO(wn, A2, S2, nu02, om2, sigma2)

# Define the log-normal function
def lognormal(x, A, a, b, c):
    return A*np.array([np.exp( (-1/(2*c**2)) * np.log((x[i] - a)/(b))**2 ) if x[i]>a else 0 for i in range(len(x))])

# sum of two log normal functions
def two_lognormals(x, A1, x01, beta1, sigma1, A2, x02, beta2, sigma2):
    return lognormal(x, A1, x01, beta1, sigma1) + lognormal(x, A2, x02, beta2, sigma2)

# voigt profile based on FWHM for Gaussian and Lorentzian components
def voigt(x, amplitude, center, fwhm_g, fwhm_l):
    sigma = fwhm_g / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation
    gamma = fwhm_l / 2  # Lorentzian HWHM
    z = ((x - center) + 1j * gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))

# sum of two Voigt profiles
def two_voigts(x, A1, t01, fwhm_g1, fwhm_l1, A2, t02, fwhm_g2, fwhm_l2):
    return voigt(x, A1, t01, fwhm_g1, fwhm_l1) + voigt(x, A2, t02, fwhm_g2, fwhm_l2)

# lorentzian function
def lorentzian(x, A, mu, FHWM):
    gamma = FHWM/2
    return A / (1 + ((x - mu) / gamma)**2)

# lorentzian function with oscillator strength
def lorentzian_area(x, f, mu, FHWM):
    gamma = FHWM/2
    # area of a lorentzian Area = A*2*Gamma
    A = f/(2*gamma)
    return A / (1 + ((x - mu) / gamma)**2)

# sum of two lorentzians
def two_lorentzians(x, A1, mu1, FWHM1, A2, mu2, FWHM2):
    return lorentzian(x, A1, mu1, FWHM1) + lorentzian(x, A2, mu2, FWHM2)

# Fcos function used to fit the coherent artifact
# more info: PhD Thesis Bastian Baudisch 2017, LMU, Riedle Group
def Fcos(t, A0, A1, tau, t0, B, phi):
    # tau here is the FWHM
    term1 = np.cos(B*(t-t0)**2 + phi)
    term2 = A0*np.exp(-4*np.log(2) * ((t - t0)**2 / (tau)**2))
    term3 = A1*((8*np.log(2))/(tau**2))*np.exp(-4*np.log(2) * ((t - t0)**2 / (tau)**2))
    return term1*(term2 - term3)

# gaussian function
def gaussian_area(t, f, t0, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    # area of a gaussian Area = A * sigma * (2*np.pi)**(0.5)
    A = f/(sigma * (2*np.pi)**(0.5))
    return A*np.exp(-(t-t0)**2/(2*sigma**2))

# gaussian function
def gaussian(t, A, t0, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    return A*np.exp(-(t-t0)**2/(2*sigma**2))

# gaussian with unit area
def gaussian_unit_area(t, t0, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    A = 1/(sigma * (2*np.pi)**(0.5))
    return A*np.exp(-(t-t0)**2/(2*sigma**2))

# sum of two gaussians
def two_gaussians(t, A1, t01, FWHM1, A2, t02, FWHM2):
    return gaussian(t, A1, t01, FWHM1) + gaussian(t, A2, t02, FWHM2)

# mono exponential decay convolved with a gaussian IRF
def gauss_conv_mono(t, A1, tau1, t0, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    term1 = (A1 / 2) * np.exp(-1 / tau1 * (t - t0))
    term2 = np.exp(0.5 * (1 / tau1**2) * sigma**2)
    term3 = 1 + erf((t - t0 - (1 / tau1) * sigma**2) / (np.sqrt(2) * sigma))
    I_N1 = term1 * term2 * term3
    return I_N1

# bi exponential decay convolved with a gaussian IRF
def gauss_conv_bi(t, A1, tau1, A2, tau2, t0, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    term1 = A1 / 2 * np.exp(-1 / tau1 * (t - t0))
    term2 = np.exp(0.5 * (1 / tau1**2) * sigma**2)
    term3 = 1 + erf((t - t0 - (1 / tau1) * sigma**2) / (np.sqrt(2) * sigma))
    
    term4 = A2 / 2 * np.exp(-1 / tau2 * (t - t0))
    term5 = np.exp(0.5 * (1 / tau2**2) * sigma**2)
    term6 = 1 + erf((t - t0 - (1 / tau2) * sigma**2) / (np.sqrt(2) * sigma))
    
    I_N2 = term1 * term2 * term3 + term4 * term5 * term6
    return I_N2

# tri exponential decay convolved with a gaussian IRF
def gauss_conv_tri(t, A1, tau1, A2, tau2, A3, tau3, t0, sigma, FWHM):
    sigma = FWHM/(2*(2*np.log(2))**0.5)
    term1 = A1 / 2 * np.exp(-1 / tau1 * (t - t0))
    term2 = np.exp(0.5 * (1 / tau1**2) * sigma**2)
    term3 = 1 + erf((t - t0 - (1 / tau1) * sigma**2) / (np.sqrt(2) * sigma))
    
    term4 = A2 / 2 * np.exp(-1 / tau2 * (t - t0))
    term5 = np.exp(0.5 * (1 / tau2**2) * sigma**2)
    term6 = 1 + erf((t - t0 - (1 / tau2) * sigma**2) / (np.sqrt(2) * sigma))
    
    term7 = A3 / 2 * np.exp(-1 / tau3 * (t - t0))
    term8 = np.exp(0.5 * (1 / tau3**2) * sigma**2)
    term9 = 1 + erf((t - t0 - (1 / tau3) * sigma**2) / (np.sqrt(2) * sigma))
    
    I_N3 = term1 * term2 * term3 + term4 * term5 * term6 + term7 * term8 * term9
    return I_N3

# sum of a gaussian and its first and second derivative
def gaussian_with_derivatives(t, A1, A2, A3, t0, FWHM):
    # Calculate sigma from FWHM
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    # Gaussian function
    gaussian = A1 * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # First derivative of Gaussian
    first_derivative = -A2 * (t - t0) / (sigma ** 2) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Second derivative of Gaussian
    second_derivative = A3 * (((t - t0) ** 2 - sigma ** 2) / (sigma ** 4)) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Sum of the Gaussian and its derivatives
    result = gaussian + first_derivative + second_derivative
    
    return result

# sum of gaussian and its second derivative 
def gaussian_and_second(t, A1, A3, t0, FWHM):
    # Calculate sigma from FWHM
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    # Gaussian function
    gaussian = A1 * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Second derivative of Gaussian
    second_derivative = A3 * (((t - t0) ** 2 - sigma ** 2) / (sigma ** 4)) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Sum of the Gaussian and its derivatives
    result = gaussian + second_derivative
    
    return result

# first derivative of a gaussian
def gaussian_first_derivative(t, A, t0, FWHM):
    # Calculate sigma from FWHM
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    return -A * (t - t0) / (sigma ** 2) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))

# second derivative of a gaussian
def gaussian_second_derivative(t, A, t0, FWHM):
    # Calculate sigma from FWHM
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    return A * (((t - t0) ** 2 - sigma ** 2) / (sigma ** 4)) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))

# sum of first and second derivative of a gaussian
def first_and_second(t, A2, A3, t0, FWHM):
    # Calculate sigma from FWHM
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    
    # First derivative of Gaussian
    first_derivative = -A2 * (t - t0) / (sigma ** 2) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Second derivative of Gaussian
    second_derivative = A3 * (((t - t0) ** 2 - sigma ** 2) / (sigma ** 4)) * np.exp(-((t - t0) ** 2) / (2 * sigma ** 2))
    
    # Sum of the Gaussian and its derivatives
    result = first_derivative + second_derivative
    
    return result

# mono exponential decay
def mono_exp(t, A, tau):
    return A*np.exp(-t/tau)

# mono exponential decay with bg
def mono_exp_with_bg(t, A, tau, bg):
    return A*np.exp(-t/tau) + bg

# bi exponential decay
def bi_exp(t, A1, tau1, A2, tau2):
    return A1*np.exp(-t/tau1) + A2*np.exp(-t/tau2)

# bi exponential decay with background
def bi_exp_with_bg(t, A1, tau1, A2, tau2, bg):
    return A1*np.exp(-t/tau1) + A2*np.exp(-t/tau2) + bg

# tri exponential decay
def tri_exp(t, A1, tau1, A2, tau2, A3, tau3):
    return A1*np.exp(-t/tau1) + A2*np.exp(-t/tau2) + A3*np.exp(-t/tau3)

# bi exponential decay with background
def tri_exp_with_bg(t, A1, tau1, A2, tau2, A3, tau3, bg):
    return A1*np.exp(-t/tau1) + A2*np.exp(-t/tau2) + A3*np.exp(-t/tau3) + bg


