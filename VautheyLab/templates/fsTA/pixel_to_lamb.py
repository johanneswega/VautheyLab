from VautheyLab.standard import *

# plot only WL spectra
debug = False

# automatic calibration
auto = True
# initial guess
p0 = [0.947, 274, 1, 1]
# set limits
lim = [0, 500]

def conversion(p, pixel, A_TA, wl_ref, A_ref):
    wl = p[0]*pixel+p[1]
    A_int = np.interp(wl_ref, wl, A_TA - p[2])
    return A_ref - p[3]*A_int

def onclick(event, x, y, ax, fig):
    if event.button == 1:  # Left mouse button click
        x.append(event.xdata)
        y.append(event.ydata)
        # Add marker at the clicked coordinates
        ax.scatter(event.xdata, event.ydata, color='k', marker='x')
        # Refresh the plot
        fig.canvas.draw()

def get_data(Ho_TA_file, WL_file, Ho_ref_file):
    WL = np.loadtxt(WL_file)    
    Ho_TA = np.loadtxt(Ho_TA_file)
    pixel = Ho_TA[:,0] 
    I_Ho = Ho_TA[:,2]
    I_WL = WL[:,2]
    # delete negative values if present
    I_Ho[I_Ho<0] = -1*I_Ho[I_Ho<0]
    I_WL[I_WL<0] = -1*I_WL[I_WL<0]
    # calculate absorbance
    A_TA = -np.log10(I_Ho/I_WL)
    # reference absorption spectrum
    Ho_ref = np.loadtxt(Ho_ref_file, skiprows=1, delimiter=',', usecols=[0,1])

    # subtract baseline and normalize
    #A_TA = (A_TA - A_TA[0])/(np.max(A_TA))

    # reference spectrum
    wl_ref = Ho_ref[:,0]
    A_ref = Ho_ref[:,1]
    #A_ref = (A_ref[wl_ref>300] - A_ref[wl_ref>300][-1]) /np.max(A_ref[wl_ref>300])
    A_ref = A_ref[wl_ref>300]
    wl_ref = wl_ref[wl_ref>300]

    return pixel, A_TA, I_Ho, I_WL, wl_ref, A_ref

def pixel_to_lambda_manual(Ho_TA_file, WL_file, Ho_ref_file):
    # get data
    pixel, A_TA, I_Ho, I_WL, wl_ref, A_ref = get_data(Ho_TA_file, WL_file, Ho_ref_file)
    
    # make first interactive figure
    fig, ax = plt.subplots(2,1,figsize=(12, 7))
    ax[0].plot(pixel, A_TA, '-r', label='Absorption TA setup')
    ax[0].set_xlabel('pixel')
    ax[0].set_ylabel('absorbance')
    ax[1].plot(wl_ref, A_ref, '-b', label='Lit. Absorption spectrum')
    ax[1].set_xlabel('wavelength / nm')
    ax[1].set_ylabel('absorbance')
    ax[0].legend()
    ax[1].legend()
    ax[0].set_title('Select peaks in TA setup spectrum!')
    x_TA = []
    y_TA = []
    cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x_TA, y_TA, ax[0], fig))
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    fig.tight_layout()
    plt.show()
    x_TA = np.array(x_TA)
    y_TA = np.array(y_TA)

    # make second interactive figure
    fig, ax = plt.subplots(2,1,figsize=(12, 7))
    ax[0].plot(pixel, A_TA, '-r', label='Absorption TA setup')
    ax[0].plot(x_TA, y_TA, 'xk')
    ax[0].set_xlabel('pixel')
    ax[0].set_ylabel('absorbance')
    ax[1].plot(wl_ref, A_ref, '-b', label='Lit. Absorption spectrum')
    ax[1].set_xlabel('wavelength / nm')
    ax[1].set_ylabel('absorbance')
    ax[0].legend()
    ax[1].legend()
    ax[0].set_title('Select the same peaks in lit. spectrum!')
    x_lit = []
    y_lit = []
    cid = fig.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x_lit, y_lit, ax[1], fig))
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    fig.tight_layout()
    plt.show()
    x_lit = np.array(x_lit)
    y_lit = np.array(y_lit)

    # do linear regression
    fig, ax = plt.subplots(1,1)
    ax.plot(x_TA, x_lit, 'ob')
    scale, shift = np.polyfit(x_TA, x_lit, 1)
    wlfine = np.linspace(0, 520, 100)
    ax.plot(wlfine, scale*wlfine + shift, '--k', label = r'px = scale $\cdot \lambda$ + shift')
    ax.set_ylabel('wavelength / nm')
    ax.set_xlabel('pixel')
    ax.set_title(r'scale = %.3g, shift = %.3g'%(scale, shift))
    fig.tight_layout()
    fig.savefig('calibration.png')
    plt.show()

    # plot absorption spectrum on the TA setup 
    fig, ax = plt.subplots(1,1)
    wl_TA = scale*pixel + shift
    ax.plot(wl_TA, A_TA/np.max(A_TA), '-r', label='TA-setup')
    ax.plot(wl_ref, A_ref/np.max(A_ref), '-b', label='literature spectrum')
    ax.legend()
    ax.set_ylabel('wavelength / nm')
    ax.set_xlabel('absorbance')
    fig.tight_layout()
    fig.savefig('comparison.png')
    plt.show()

def pixel_to_lambda(Ho_TA_file, WL_file, Ho_ref_file):
    # get data
    pixel, A_TA, I_Ho, I_WL, wl_ref, A_ref = get_data(Ho_TA_file, WL_file, Ho_ref_file)

    # make plot
    fig, ax = plt.subplots(ncols=1,nrows=1)
    if debug==True:
        plt.plot(pixel, I_Ho, '-r')
        plt.plot(pixel, I_WL, '-b')
        plt.show()

    # least squares fit
    res = least_squares(conversion, x0=p0,
                        args=(pixel[(pixel<lim[1])&(pixel>lim[0])], 
                              A_TA[(pixel<lim[1])&(pixel>lim[0])], wl_ref, A_ref))
    p = res.x
    scale = p[0]
    shift = p[1]

    # plot results
    print('scale = %.3g, shift = %.3g'%(scale,shift))

    # write results to file 
    wl = scale*pixel+shift
    # Open the file in write mode
    with open('wl.txt', 'w') as file:
        # Write each value in the list as a separate line in the file
        for item in wl:
            file.write(str(item) + '\n')

    ax.plot(scale*pixel+shift, A_TA, '-b', linewidth=1.5, label='TA setup')
    ax.plot(wl_ref, A_ref, '-r', linewidth=1.5, label='Ref.')
    ax.set_xlabel(r'$\lambda /$ nm')
    ax.set_ylabel(r'norm. Abs.')
    ax.set_xticks(np.linspace(300,800,6))
    ax.set_yticks([])
    ax.legend()
    ax2 = ax.secondary_xaxis("top", functions=(lambda x: (x-shift)/scale,lambda x: x*scale + shift))
    ax2.set_xlabel('pixel')
    ax2.set_xticks(np.linspace(0,500,6))
    fig.tight_layout()
    fig.savefig('pixel_to_lambda.png')

# get data files
files = [f for f in os.listdir() if f.endswith('.dat')]
for i in range(len(files)):
    if 'WL' in files[i]:
        WL = files[i]
    if 'HOLMIUM' in files[i]:
        HO = files[i]

if auto==True:
    pixel_to_lambda(HO, WL, 'HOLM_ref.csv')
else:
    pixel_to_lambda_manual(HO, WL, 'HOLM_ref.csv')
plt.show()