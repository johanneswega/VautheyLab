from VautheyLab.miscellaneous import moving_average as ma
from VautheyLab.transient_absorption import *
from VautheyLab.standard import *
from VautheyLab.standard_figures import *

def plot_data(file, label, color, scat):
    data = np.loadtxt(file, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    if scat!=None:
        y[(x >= scat[0]) & (x <= scat[1])] = np.nan
    x = (1/x)*10**4
    y = -y/np.min(y)
    xm = ma(x, 5)
    ym = ma(y, 5)
    #norm = -1*np.nanmin(ym[xm>24])
    ax.plot(x, y, '-', alpha=0.1)
    ax.plot(xm, ym, '-', label=label, color=color)


specs = ['A', 'B', 'C', 'D']
solvs = ['0 mM', '1000 mM']
colors = ['b', 'r']
scatter = [None, None]

for k in range(len(specs)):
    fig, ax  = figure_spectum()
    spec = specs[k]

    for i in range(len(solvs)):
        data = np.loadtxt('%s/lifetimes.txt'%(solvs[i]), delimiter=',')
        taus = data[:,0]
        tau = taus[specs.index(spec)]
        if tau<1e3:
            label = '%s (%.3g ns)'%(solvs[i], tau)
        else:
            label = '%s (%.3g µs)'%(solvs[i], tau/1e3)
        plot_data('%s/%s.txt'%(solvs[i], spec), label, colors[i], scatter[i])

    ax.set_ylabel(r'norm. $\Delta A$')
    ax.legend()
    title = [r'$\textbf{%s} \rightarrow $'%s if s == spec else r'$\text{%s} \rightarrow $'%s for s in specs]
    ax.set_title(r''.join(title))
    fig.tight_layout()
    plt.show()