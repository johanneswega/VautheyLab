from VautheyLab.standard import *
from VautheyLab.miscellaneous import moving_average as ma, find_index

# function to get data 
def get_data(file, f_type):
    if f_type=='r_com':
        data = np.loadtxt(file, skiprows=24)
        t = data[:,0]
        r_com = data[:,1]
        return t, r_com

# function to plot / compare trajectories
def plot_trajectory(files, colors, labels, f_type, figsize=(6, 3.5), find_time=None,
                    xlim=False, MA=False, MA_npoints=10, outside=True, linewidth=1, save=None):
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if MA==False:
            ax.plot(x, y, '-', color=colors[i], label=labels[i], linewidth=linewidth)
        else:
            ax.plot(x, y, '-', color=colors[i], alpha=0.2, linewidth=linewidth)
            ax.plot(ma(x,MA_npoints), ma(y,MA_npoints), '-', color=colors[i], label=labels[i], linewidth=linewidth)
    # find times for certain width
    if find_time!=None:
        for i in range(len(find_time)):
            print(find_time[i], round(x[find_index(y, find_time[i])]))
            ax.plot(x[find_index(y, find_time[i])], y[find_index(y, find_time[i])], 'or', markersize=5)
    if f_type=='r_com':
        ax.set_ylabel(r'$r_{\text{com}} / $ nm')
        ax.set_xlabel(r'$t / $ ps')
    if xlim!=None:
        ax.set_xlim(xlim)
    if outside==True:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    else:
        ax.legend()
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# function to plot / compare trajectories using a stack
def plot_trajectory_stack(files, colors, labels, f_type, figsize=None, save=None, ylim=False,
                          xlim=False, MA=False, MA_npoints=10, outside=True, linewidth=1):
    # make figure
    if figsize==None:
        figsize = (6, len(files)*3)
    fig, ax = plt.subplots(len(files), 1, figsize=figsize, sharex=True, sharey=True)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        if MA==False:
            ax[i].plot(x, y, '-', color=colors[i], label=labels[i], linewidth=linewidth)
        else:
            ax[i].plot(x, y, '-', color=colors[i], alpha=0.2, linewidth=linewidth)
            ax[i].plot(ma(x,MA_npoints), ma(y,MA_npoints), '-', color=colors[i], label=labels[i], linewidth=linewidth)
        if f_type=='r_com':
            ax[i].set_ylabel(r'$r_{\text{com}} / $ nm')
        if outside==True:
            ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        else:
            ax[i].legend()
    ax[-1].set_xlabel(r'$t / $ ps')
    if xlim!=None:
        ax[-1].set_xlim(xlim)
    if ylim!=False:
        ax[-1].set_ylim(ylim)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# histogram plotter
def histogram(files, colors, labels, f_type, nbins, width, figsize=(6, 3.5), outside=True, ylim=False, save=None):
    # make figure
    fig, ax = plt.subplots(1,1,figsize=figsize)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        ax.hist(y, nbins, color=colors[i], width=width, edgecolor='k', label=labels[i])
    if f_type=='r_com':
        ax.set_xlabel(r'$r_{\text{com}} / $ nm')
    ax.set_ylabel(r'$N$')
    if outside==True:
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    else:
        ax.legend()
    if ylim!=False:
        ax.set_ylim(ylim)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()

# histogram plotter
def histogram_stack(files, colors, labels, f_type, nbins, width, figsize=None, outside=True, ylim=False, save=None):
    # make figure
    if figsize==None:
        figsize = (6, len(files)*3)
    fig, ax = plt.subplots(len(files), 1, figsize=figsize, sharex=True, sharey=True)
    # get data and plot
    for i in range(len(files)):
        x, y = get_data(files[i], f_type)
        ax[i].hist(y, nbins, color=colors[i], width=width, edgecolor='k', label=labels[i])
        ax[i].set_ylabel(r'$N$')
        if outside==True:
            ax[i].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        else:
            ax[i].legend()
    if f_type=='r_com':
        ax[-1].set_xlabel(r'$r_{\text{com}} / $ nm')
    if ylim!=False:
        ax[-1].set_ylim(ylim)
    fig.tight_layout()
    if save!=None:
        fig.savefig(save, transparent=True)
    plt.show()