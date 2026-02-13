# VautheyLab

**Object-oriented Python toolkit for spectroscopic data analysis** for steady-state and time-resolved spectroscopy, cyclic voltammetry, and QM/MD workflows developed during my PhD at the Vauthey Group (University of Geneva).  

---

## ⚙️ Overview

This package provides modular, object-oriented tools for handling and analyzing experimental data. Objects are created for each type of experiment or analysis technique, while the mathematical processing and analysis happens in the background.  

Templates for common data analysis workflows are available in `VautheyLab/templates`, and some example scripts can be found in the `examples` folder.  

> Note: Documentation is currently sparse. For more detailed guidance on the capabilities of the package, best to talk to me directly. I can show you how to use it and explain what happens under the hood.

---

## 📦 Installation

It is recommended to use **Python 3.10** in a virtual environment (e.g., via `venv` or `conda`).

1. Clone the repository:

```bash
git clone https://github.com/johanneswega/VautheyLab.git
cd VautheyLab
```

2. Install the package in editable mode:

```bash
pip install -e .
```

The -e flag (--editable) allows you to modify the scripts locally while using the package. As I use LaTeX for nice mathematical axes labels it is also necessary to have LaTeX installed on your machine.

---

## 🛠 Usage Example

For example, if you measured an absorption spectrum on our Cary50 spectrometer (`abs_file.csv`), you can plot it using the `Absorption` class from `VautheyLab.steady_state`:

```python
from VautheyLab.steady_state import Absorption

a = Absorption(
    files=['abs_file.csv'],
    cuts=[(200, 800)],
    units='wn',
    colors=['r'],
    norm=True,
    labels=['my molecule in X']
)
a.show()
```
In almost all class objects, arguments like `files`, `cuts`, `colors`, and `labels` are given as **lists**. This makes it easy to plot and compare multiple spectra in the same figure.

You can also directly access the loaded data and analyze it with your own custom code. For example, the wavelength, wavenumber, and absorbance of the *i-th* file are stored as NumPy arrays:
```python
wavelength = a.wl[0]
wavenumber = a.wn[0]
absorbance = a.A[0]
```

However, most of the classes are not just for plotting. Most of them come with built-in analysis methods. Say you want to estimate the oscillator strength and radiative rate constant for the $S_1 \leftarrow S_0$ transition by integrating the absorption spectrum between 440–600 nm using the Strickler-Berg analysis. For this, xou just need to pass two extra arguments when creating the Absorption object:
- `conc` → concentration in M
- `pathlength` → optical pathlength in cm

With these two additional init arguments, the class will directly calculate the extinction spectrum. To calcualte the oscillator strength you can then call `calc_oscillator_strength()` method of the class: 

```python
from VautheyLab.steady_state import Absorption

# concentration in M
c = 24e-6
# pathlength in cm
l = 1

# initialize absorption class
a = Absorption(files=['abs_file.csv'],
             cuts=[(350, 800)],
             colors=['r'],
             units='wn',
             conc=[24e-6], 
             pathlegnths=[1],
             labels=['your molecule'])

# use calc_oscillator_strength(limits, n, nu0, file_index)
# integration limits as list in nm
# refractive index 
# center frequency in kk
a.calc_oscillator_strength([440, 600], 1.421, 18.8, 0)
a.show()
```

This will automatically calculate and output the following useful properties: 

- The oscillator strenth calculated over the region is: $f = 0.13$
- The TDM is: $\mu_{\text{TDM}} = 3.75\,\text{D}$
- Radiative rate constant $k_{\text{rad}} = 5.81 \times 10^8 \,\text{s}^{-1}$
- Radiative lifetime $\tau_\text{rad} = 17.2\,\text{ns}$

Of course this is not the only method of the `absorption` class and a plethora of other methods are implemented like:
- `get_concentration` (get concentration of sample providing exctinction coefficeint)
- `find_max` (find wavelength of maximum absorbance)
- `find` (find abs at a certain wavelength) 
- `plot_calculated` (compare with Gaussian output convolved calculated spectrum / btw the class works also for FT-IR spectra)
- `solvchrom` (make a solvatochromic plot of all files when solvent names are given as attributes, solvent parameters are directly extracted)
- `plot_diff` (plot difference spectra)
...

So the modular approch is pretty versatile I would say. I don’t have time (yet) to fully document every module. As you can see this section alone is already long just for absorption. If you have questions and problems either ask me or study how things are done under the hood by looking at the implementations of the respective class objects. For instance, the `absorption`class is implemented in `VautheyLab/absorption.py`. In this file you’ll find:

- all initialization arguments
- all analysis methods of the class

and can see how things are implemented under the hood. Comments are provided in the code.

Whenever I get time, I’ll keep adding more examples for different experiments/analysis routines in the `VautheyLab/examples` folder. 