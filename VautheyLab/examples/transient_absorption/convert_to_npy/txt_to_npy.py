from VautheyLab.transient_absorption import txt_to_npy

# will export .npy file as dA_from_txt.npy
txt_to_npy(t_file='time.txt', wl_file='wavelength.txt', dA_file='TA.txt')
