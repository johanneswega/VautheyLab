# VautheyLab

**Object-oriented Python toolkit for spectroscopic data analysis** — steady-state and time-resolved spectroscopy, cyclic voltammetry, and QM/MD workflows — developed during my PhD at the Vauthey Group, University of Geneva.  

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

The -e flag (--editable) allows you to modify the scripts locally while using the package.

---

🛠 Usage Example

For example, if you measured an absorption spectrum on the Cary50 spectrometer (`abs_file.csv`), you can plot it using the `Absorption` class from `VautheyLab.steady_state`:

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