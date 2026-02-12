from setuptools import setup, find_packages

setup(
    name='VautheyLab',
    version='0.1',
    packages=find_packages(),  # Automatically finds your 'VautheyLab' folder and submodules
)

# to install it in editable mode 
# pip3 install -e .

# to install it permanently so changed do not apply if you modify a file 
# pip3 install .