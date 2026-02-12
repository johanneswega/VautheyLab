from VautheyLab.standard import *
from VautheyLab.miscellaneous import solvs, find_index

class Compare_traces:
    def __init__(self, files, irf_file, figsize=None, outside=False, fontsize=10, plotfit=False, taus=None,
                 ylim=None, xlim=None, xticks=True, cleanup=5, notime=False, irfnotime=False, ulf=False,
                 yticks=True, colors=['r'], labels=[''], tightlayout=True, savefig=False, yscale='log', waterfall=0):
        
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.files = files
        # decide whether to normalize or not
        self.irf_file = irf_file
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        # specify list of colors
        self.colors = colors
        self.ulf = ulf
        # specify list of labels
        self.labels = labels
        self.cleanup = cleanup
        self.taus = taus
        # tight layout
        self.tightlayout = tightlayout
        self.savefig = savefig
        self.plotfit = plotfit
        # yscale
        self.yscale = yscale
        self.outside = outside
        self.fontsize = fontsize
        self.notime = notime
        self.irfnotime = irfnotime
        self.waterfall = waterfall
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # plot data
        self.plot()

    def read_data(self, fname, notime):
        if self.ulf == False:
            if notime==False:
                fh = open(fname, 'r')
                counts = []
                t = []
                for line in fh: 
                    if line[0]!=' ' and line[0]!='*':
                        line = line.split()
                        if line[1]=='0':
                            line[1] = line[1].replace('0','1')
                        t.append(float(line[0]))
                        counts.append(float(line[1]))
            else:
                fh = open(fname, 'r')
                counts = []
                t = []
                for line in fh: 
                    if line[0]!=' ' and line[0]!='*':
                        line = line.replace('0','1')
                        counts.append(float(line))
                # for normal length = corresponds to 50 ns in total 
                if len(counts)==4096:
                    dt = 50/(len(counts) -1)
                    t = [i*dt for i in range(len(counts))]
            return np.array(t), np.array(counts)
        else:
            data = np.loadtxt(fname, skiprows=10)
            I = data
            # ns per channel 
            t = np.array([i*0.0040 for i in range(len(I))])  
            return t, I
    
    def read_data_fit(self, fname):
        data = np.loadtxt(fname, skiprows=4, delimiter=',')
        return data[:,0], data[:,1], data[:,2]
    
    # function to extract taus from export fit file
    def get_tau(self):
        taus = []
        for i in range(len(self.files)):
            fh = open(self.files[i], 'r')
            for line in fh:
                if line.startswith('#') and 'lifetime' in line:
                    tau = float(line.split()[4][1:])
                    taus.append(tau)
            fh.close()
        return taus
    
    def comp_taus(self, ylim=None, solv=None):
        # comparision of all taus
        fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
        ticks = [i for i in range(len(self.colors))]
        for i in range(len(self.colors)):
            ax.plot(i, self.taus[i], 's', color=self.colors[i])

        ax.set_xticklabels(self.labels)
        ax.set_xticks(ticks)
        if ylim!=None:
            ax.set_ylim(ylim)
        ax.set_ylabel(r'$\tau / \text{ns}$')
        fig.tight_layout()
        fig.savefig('taus.svg', transparent=True)

        # against df and eta 
        if solv!=None:
            # get viscosity
            eta = [solvs[i][3] for i in solv]
            # get df
            df = [solvs[i][2] for i in solv]

            # df vs tau
            fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
            for i in range(len(self.colors)):
                ax.plot(df[i], self.taus[i], 's', color=self.colors[i], label=solv[i])
            ax.set_ylabel(r'$\tau / \text{ns}$')
            ax.set_xlabel(r'$\Delta f$')
            ax.legend()
            fig.tight_layout()
            fig.savefig('taus_vs_df.svg', transparent=True)

            # df vs tau
            fig, ax = plt.subplots(1,1,figsize=(5, 3.5))
            for i in range(len(self.colors)):
                ax.plot(eta[i], self.taus[i], 's', color=self.colors[i], label=solv[i])
            ax.set_ylabel(r'$\tau / \text{ns}$')
            ax.set_xlabel(r'$\eta$ / cP')
            ax.set_xscale('log')
            ax.legend()
            fig.tight_layout()
            fig.savefig('taus_vs_eta.svg', transparent=True)
            
    def plot(self):
        tIRF, IRF = self.read_data(self.irf_file, self.irfnotime)
        IRF[IRF<self.cleanup] = 0

        for i in range(len(self.files)):
            t, counts, fit = self.read_data_fit(self.files[i])
            # reject all data with less than 1 counts
            counts[counts==0] = 1
            # redefine time zero to be on the peak of the IRF
            t = t - t[counts==np.max(counts)][0]
            # get background
            bg = np.mean(counts[t<-1])
            counts = counts - bg
            norm = np.max(fit)
            counts = counts/norm
            IRF = IRF/np.max(IRF)
            shift = t[fit==np.max(fit)]
            t0 = t[fit==np.max(fit)]
            shift = t[t<t0][find_index(fit[t<t0], np.max(fit)/2)]
            t = t - shift
            if self.plotfit==True:
                fit = fit - bg
                fit = fit/norm
                self.ax.plot(t, counts + i*self.waterfall, '.', color=self.colors[i], alpha=0.1)
                self.ax.plot(t, fit + i*self.waterfall, '-', color=self.colors[i],  label=self.labels[i])
                #self.ax.fill_between(tIRF - tIRF[IRF==np.max(IRF)], i*self.waterfall, IRF + i*self.waterfall, color='gray', alpha=0.1, label='IRF')
            else:
                #self.ax.plot(t, counts, '-', color=self.colors[i], label=self.labels[i])
                self.ax.plot(t, counts, '.', color=self.colors[i], alpha=0.1)
        self.ax.fill_between(tIRF - tIRF[IRF==np.max(IRF)], 0, IRF, color='gray', alpha=0.1, label='IRF')

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)

        self.ax.set_ylabel(r'norm. counts')
        self.ax.set_xlabel(r'$\Delta t$ / ns')
        self.ax.legend()
        self.ax.set_yscale(self.yscale)

        if self.outside==True:
            self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
        else:
            self.ax.legend(fontsize=self.fontsize)
        if self.yticks == False:
            self.ax.set_yticks([])
        if self.xticks == False:
            self.ax.set_xticks([])
        if self.ylim != None:
            self.ax.set_ylim(self.ylim)
        if self.xlim != None: 
            self.ax.set_xlim(self.xlim)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            self.fig.savefig('trace_comp.svg', transparent=True)
        plt.show()
