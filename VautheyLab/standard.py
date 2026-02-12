import numpy as np
from matplotlib import pyplot as plt
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))
from matplotlib.animation import FuncAnimation, writers
from scipy.optimize import curve_fit, least_squares
import scipy.constants as sc