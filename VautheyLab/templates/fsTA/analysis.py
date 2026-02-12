from VautheyLab.standard import *
import os
import shutil
from tqdm import tqdm
from PIL import Image
from PyPDF2 import PdfMerger

# define excitation wavelength in nm
ex_wl = 400
# define scans to be excluded
exclude = None
# define scale
vmin = -20
vmax = +20
# show
show = True
# chirp correction with or without OKE
ok = False

def figshow(figure):
    for i in plt.get_fignums():
        if figure != plt.figure(i):
            plt.close(plt.figure(i))
    plt.show()

def fit(l, a, b, c):
    return a + 10**5 * b * l**-2 + 10**6 * c * l**-4

def get_individual_dA(TA_file):
    fh = open(TA_file)
    number_scans = 0
    timedelays = 0
    skip = 0
    for line in fh:
        if 'scans' in line:
            number_scans += int(line[line.find('=')+2])
        if 'time-delays' in line: 
            timedelays += int(line[line.find('=')+2:])
        if '%' in line:
            skip += 1
    fh.close()
    scans = []
    for i in range(number_scans):
        scans.append([])
    df = np.loadtxt(TA_file, skiprows=skip)
    # time in ps
    t = df[:,0][:timedelays]*10**12
    df = np.vsplit(df, number_scans)
    # get indices of error columns
    index = [0, 1]
    for i in range(len(df[0][0,:])-2):
        if not i%2==0:
            index.append(i+2)
    # delete theese columns
    for i in range(number_scans):
        df[i] = np.delete(df[i], index, axis=1)
        # convert from micro OD to mOD
        df[i] = df[i]*10**-3
    return t, df

def plot_2D_map(t, l, dA, ax, sinh=None):
    l = (1/l)*10**4
    if sinh==None:
        ax.contour(l, t, dA, levels=15, linewidths=0.12, colors='k')
        D = ax.pcolormesh(l, t, dA, vmin=vmin, vmax=vmax, cmap="RdBu_r")
        cbar = plt.colorbar(D, ax=ax)
        cbar.set_label(r'$\Delta{A} \, \times 10^{3}$')  
    else:
        ax.contour(l, t, np.arcsinh(dA), levels=15, linewidths=0.12, colors='k')
        D = ax.pcolormesh(l, t, np.arcsinh(dA), vmin=np.arcsinh(vmin), vmax=np.arcsinh(vmax), cmap="RdBu_r")
        cbar = plt.colorbar(D, ax=ax)
        cbar.set_label(r'arcsinh($\Delta{A} \, \times 10^{3}$)')    
    ax.set_xlabel(r'$\tilde{\nu} \, \times 10^{-3} \, cm^{-1}$')
    ax.set_ylabel(r'$\Delta t /$ ps')
    ax.set_yscale('symlog')
    ax.set_yticks([0, 1, 2, 10, 100, 1000])
    ax.set_yticklabels([0, 1, 2, 10, 100, 1000])
    ax.set_ylim([-0.5, t[-1]])
    ax.invert_xaxis()
    f = lambda x: (1/x)*10**-4
    g = lambda x: (1/x)*10**+4
    ax2 = ax.secondary_xaxis("top", functions=(g,f))
    ax2.set_xlabel(r'$\lambda /$ nm')
    
def plot_background_spectra(ax, l, dA, i):
    l = (1/l)*10**4
    ax.plot(l, dA, '-', label=r'scan %i'%(i+1))
    ax.set_ylabel(r'$\Delta{A} \, \times 10^{3}$')  
    ax.set_xlabel(r'$\tilde{\nu} \, \times 10^{-3} \,$ cm$^{-1}$')
    ax.invert_xaxis()
    
def plot_overview(ax, t, l, dA, delay):
    l = (1/l)*10**4
    c = plt.cm.gist_rainbow(np.linspace(0, 1, len(delay)))
    for i in range(len(delay)):
        index = min(range(len(t)), key=lambda j: abs(t[j]-delay[i]))
        if delay[i]>=1000:
            lab='%.3g ns'%(delay[i]/1000)
        else:
            lab='%.3g ps'%(delay[i])
        ax.plot(l, dA[index,:], '-', linewidth=1, color=c[i], label=lab)
    ax.invert_xaxis()
    ax.set_xlim([29, 13])
    ax.set_xticks(np.linspace(28,14,8))
    index = min(range(len(t)), key=lambda j: abs(t[j]-3))
    ax.set_ylim([vmin, vmax])
    ax.axhline(y=0, color='k', linewidth=1)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
    ax.set_ylabel(r'$\Delta{A} \, \times 10^{3}$')  
    ax.set_xlabel(r'$\tilde{\nu} \, \times 10^{-3} \,$ cm$^{-1}$')
    f = lambda x: (1/x)*10**-4
    g = lambda x: (1/x)*10**+4
    ax2 = ax.secondary_xaxis("top", functions=(g,f))
    ax2.set_xlabel(r'$\lambda$ / nm')

def onclick(event, x, y, ax, fig):
    if event.button == 1:  # Left mouse button click
        x.append(event.xdata)
        y.append(event.ydata)
        # Add marker at the clicked coordinates
        ax.scatter(event.xdata, event.ydata, color='red', marker='o')
        # Refresh the plot
        fig.canvas.draw()

def load_TA_data(ex_wl, TA_file, exclude):
    # progress bar
    progress_bar = tqdm(total=7)

    # make results folder 
    if os.path.exists('Results')==True:
        shutil.rmtree('Results')
    os.mkdir('Results')

    # get wavelength axis
    wl = np.loadtxt('wl.txt') 
    # load dataframe
    t, df = get_individual_dA(TA_file)

    # plot individual scans
    nscans = len(df)
    fig1,ax1 = plt.subplots(ncols=1,nrows=nscans)
    fig1.set_size_inches(6, 4*nscans)
    for i in range(nscans):
        plot_2D_map(t, wl, df[i], ax1[i])
        ax1[i].set_title(r'scan %i'%(i+1))
    fig1.tight_layout()
    fig1.savefig("Results/01_raw_scans.png")
    progress_bar.update(1)

    # plot background spectra
    fig2, ax2 = plt.subplots(ncols=1, nrows=1)
    fig2.set_size_inches(6, 4)
    for i in range(nscans):
        bg = np.mean(df[i][t<-1.8,:], axis=0)
        # subtract background
        for j in range(len(t)):
            df[i][j,:] = df[i][j,:] - bg 
        plot_background_spectra(ax2, wl, bg, i)
    ax2.legend()
    ax2.axvline(x=(1/ex_wl)*10**4, color='k')
    f = lambda x: (1/x)*10**-4
    g = lambda x: (1/x)*10**+4
    ax22 = ax2.secondary_xaxis("top", functions=(g,f))
    ax22.set_xlabel(r'$\lambda$ / nm')
    fig2.tight_layout()
    fig2.savefig("Results/02_background.png")
    progress_bar.update(1)

    # plot background corrected spectra
    fig3,ax3 = plt.subplots(ncols=1,nrows=nscans)
    fig3.set_size_inches(6, 4*nscans)
    for i in range(nscans):
        plot_2D_map(t, wl, df[i], ax3[i])
        ax3[i].set_title(r'scan %i - bg corrected'%(i+1))
    fig3.tight_layout()
    fig3.savefig("Results/03_raw_scans_bg_corrected.png")  
    progress_bar.update(1)

    # calculate and plot averaged spectrum
    # exclude traces if wanted
    if exclude!=None:
        df_new = []
        for i in range(nscans):
            if not i in exclude:
                df_new.append(df[i])
        dA = np.mean(df_new, axis=0)
    else:
        dA = np.mean(df, axis=0)
    fig4,ax4 = plt.subplots(ncols=1, nrows=1)
    fig4.set_size_inches(6, 4)
    plot_2D_map(t, wl, dA, ax4)
    fig4.tight_layout()
    fig4.savefig("Results/04_averaged_bg_corrected.png")
    progress_bar.update(1)

    # chirp correct data 
    # plot chirp line
    fig5,ax5 = plt.subplots(ncols=1, nrows=1)
    # retrieve OKE data
    if ok==True:
        fig5.set_size_inches(6, 4)
        p = np.loadtxt('ok.txt')
        plot_2D_map(t, wl, dA, ax5)
        ax5.plot(np.linspace(wl[0],wl[-1], 1000)**-1 *10**4,fit(np.linspace(wl[0],wl[-1], 1000), *p), '--k')
        fig5.tight_layout()
        fig5.savefig("Results/05_Data_chirp_line.png")
        progress_bar.update(1)
    # get fit line by hand
    else:
        fig5.set_size_inches(8, 8)
        plot_2D_map(t, wl, dA, ax5)
        x = []
        y = []
        cid = fig5.canvas.mpl_connect('button_press_event', lambda event: onclick(event, x, y, ax5, fig5))
        figshow(fig5)
        x = np.array(x)
        y = np.array(y)
        p, pcov = curve_fit(fit, (1/x)*10**4, y)

        fig51,ax51 = plt.subplots(ncols=1, nrows=1)
        plot_2D_map(t, wl, dA, ax51, False)
        ax51.plot(np.linspace(wl[0],wl[-1], 1000)**-1 *10**4,fit(np.linspace(wl[0],wl[-1], 1000), *p), '--k')
        fig51.tight_layout()
        fig51.savefig("Results/05_Data_chirp_line.png")
        progress_bar.update(1)
        figshow(fig51)

    # plot chirp corrected 2D TA
    for i in range(len(wl)):
        t0 = fit(wl[i], *p)
        dA[:,i] = np.interp(t, t-t0, dA[:,i])            
    fig6,ax6 = plt.subplots(ncols=2, nrows=1)
    fig6.set_size_inches(10, 4)
    plot_2D_map(t, wl, dA, ax6[0])
    plot_2D_map(t, wl, dA, ax6[1], False)
    ax6[0].axhline(y=0, color='k')
    ax6[1].axhline(y=0, color='k')
    fig6.tight_layout()
    fig6.savefig("Results/06_Data_chirp_corrected.png")
    progress_bar.update(1)

    # plot overview spectrum
    fig7,ax7 = plt.subplots(ncols=1, nrows=1)
    fig7.set_size_inches(5, 3.5)
    delay = [-2, -1, 0.5, 1, 2, 5, 10, 50, 100, 250, 500, 1000, 1200, 1800]
    plot_overview(ax7, t, wl, dA, delay)
    fig7.tight_layout()
    fig7.savefig("Results/07_overview.png")
    progress_bar.update(1)
    
    # generate results file
    if os.path.exists('Results.pdf')==True:
            os.remove('Results.pdf')
    names = os.listdir('Results')
    names = np.sort(names)
    for i in range(len(names)):
        image = Image.open('Results/%s'%names[i])
        im = image.convert('RGB')
        im.save('Results/%s.pdf'%(names[i][:-4]))
    merger = PdfMerger()
    for i in range(len(names)):
        merger.append('Results/%s.pdf'%(names[i][:-4]))
    merger.write('Results.pdf') 
    merger.close()
    np.save('dA.npy', np.array([t, wl, dA], dtype=object))
    pdat = np.zeros((len(t)+1, len(wl)+1))
    pdat[1:,0] = t
    pdat[0, 1:] = wl
    pdat[1:, 1:] = dA
    np.savetxt('dA.pdat', pdat, header='ps*nm', delimiter=',')

    progress_bar.close()

    if show==True:
        plt.show()

# get data files
files = [f for f in os.listdir() if f.endswith('.dat')]
for i in range(len(files)):
    if 'SAMPLE' in files[i]:
        sample = files[i]

load_TA_data(
             ex_wl,
             sample,
             exclude)
