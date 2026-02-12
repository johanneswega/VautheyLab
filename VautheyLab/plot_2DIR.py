from VautheyLab.standard import *
from VautheyLab.miscellaneous import find_index

def round_sig(n):
    if n == 0:
        return 0  # Edge case: 0 remains 0
    power = np.floor(np.log10(abs(n)))  # Find the power of 10
    factor = 10 ** power  # Compute the rounding factor
    return np.ceil(n / factor) * factor  # Round up and scale back

class twoDIR:
    def __init__(self, file, delay, pump_cuts=None, probe_cuts=None, norm=False, normat=None, minnorm=False, ylim=None, xlim=None, xticks=True,
                 yticks=True, export=False, tightlayout=True, cmap='seismic', lines=True, diagonal=True, flip=False,
                 savefig=False, figsize=None, scale=None, nlevels=50, white=True):          
        # set relevant parameters
        self.figsize = figsize
        # get files
        self.file = file
        # get wavelength cuts
        self.delay = delay
        self.pump_cuts = pump_cuts
        self.probe_cuts = probe_cuts
        # decide whether to normalize or not
        self.norm = norm
        # if you want to norm on GSB
        self.minnorm = minnorm
        # to normalize at a specific wl
        self.normat = normat
        if self.normat!=None and self.normat[1]<100:
            self.normat[1] = (1/self.normat[1])*10**4
        # specify x and y limits for plot
        self.ylim = ylim
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        # whether to export data or not
        self.export = export
        # contour specific things
        self.nlevels = nlevels
        self.white = white
        self.scale = scale
        self.cmap = cmap
        self.lines = lines
        self.diagonal = diagonal
        self.flip = flip
        self.tightlayout = tightlayout
        # save figure
        self.savefig = savefig
        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data()
        self.plot()

    def read_data(self):
        # get t2 delays
        with open(self.file, "r") as file:
            t2 = file.readline().split(',')
            t2 = [float(x) for x in t2]
        self.t2 = np.array(t2[2:])
        # print t2 delays
        print("\n The 2D-Data set contains the following t2-delays: \n")
        for i in range(len(self.t2)):
            print("t2 = %.3g ps \n"%(self.t2[i]))

        data_original = np.loadtxt(self.file, skiprows=1, delimiter=',')
        # split into equal chunks of 64 pixels
        data = np.split(data_original, indices_or_sections=np.shape(data_original)[0] // 64, axis=0)  
        # make empty data frame
        self.dA2D = np.zeros((np.shape(data)[0], np.shape(data)[-1], np.shape(data)[1]))
        # extract probe axis
        self.nu3 = data[0][:,0]
        # extract pump axis 
        self.nu1 = []
        for i in range(np.shape(data)[0]):
            self.nu1.append(data[i][0,1])

        # fill in data frame
        for i in range(len(data_original[:,0])):
            # find pump and probe index value
            probe_index = find_index(self.nu3, data_original[i,:][0])
            pump_index = find_index(self.nu1, data_original[i,:][1]) 
            # go through row and put in data frame
            row = data_original[i,:][2:]
            for j in range(len(self.t2)):
                self.dA2D[pump_index, j, probe_index] = row[j]

        # cut data if wanted
        self.nu1 = np.array(self.nu1)
        self.nu3 = np.array(self.nu3)
        if self.probe_cuts!=None:
            dA2D_cut = np.zeros((len(self.nu1), len(self.t2), len(self.nu3[(self.nu3>self.probe_cuts[0])&(self.nu3<self.probe_cuts[1])])))
            for i in range(len(self.t2)):
                dA2D_cut[:, i, :] = self.dA2D[:, i, (self.nu3>self.probe_cuts[0])&(self.nu3<self.probe_cuts[1])]
            self.nu3 = self.nu3[(self.nu3>self.probe_cuts[0])&(self.nu3<self.probe_cuts[1])]
            self.dA2D = dA2D_cut

        if self.pump_cuts!=None:
            dA2D_cut = np.zeros((len(self.nu1[(self.nu1>self.pump_cuts[0])&(self.nu1<self.pump_cuts[1])]), len(self.t2), len(self.nu3)))
            for i in range(len(self.t2)):
                dA2D_cut[:, i, :] = self.dA2D[(self.nu1>self.pump_cuts[0])&(self.nu1<self.pump_cuts[1]), i, :]
            self.nu1 = self.nu1[(self.nu1>self.pump_cuts[0])&(self.nu1<self.pump_cuts[1])]
            self.dA2D = dA2D_cut

    def get_anharmonicity(self, range_pump, range_probe):
        dAlook = self.dA2D[:, find_index(self.t2, self.delay), :]  
        dAlook = dAlook[(self.nu1>range_pump[0])&(self.nu1<range_pump[1]), :]
        nu1now = self.nu1[(self.nu1>range_pump[0])&(self.nu1<range_pump[1])]    
        dAlook = dAlook[:, (self.nu3>range_probe[0])&(self.nu3<range_probe[1])]
        nu3now = self.nu3[(self.nu3>range_probe[0])&(self.nu3<range_probe[1])]  
        max_y, max_x = np.where(dAlook == np.max(dAlook))
        print("Max indices:", max_y, max_x)
        print(nu3now[max_x], nu1now[max_y])
        self.ax.plot(nu3now[max_x], nu1now[max_y], 'xr', markersize=15)
        self.ax.plot(nu3now[max_x], nu1now[max_y], 'or', markersize=5)
        self.ax.contour(nu3now, nu1now, dAlook, colors='k', levels=10)
        

    def plot(self):
        # find desired t2 index
        self.dA = self.dA2D[:, find_index(self.t2, self.delay), :]
        # fix scale
        if self.scale==None:
            self.scale = [-1*round_sig(np.max(self.dA)), round_sig(np.max(self.dA))]
        # levels and contours
        levels = np.linspace(self.scale[0], self.scale[1], self.nlevels)
        colors = plt.cm.seismic(np.linspace(0, 1.0, self.nlevels))
        if self.white==True:
            colors[int(self.nlevels/2)] = 0
            colors[int(self.nlevels/2)-1] = 0
            colors[int(self.nlevels/2)-2] = 0

        print(np.shape(self.dA), np.shape(self.nu1), np.shape(self.nu3))

        # plot
        if self.flip==False:
            D = self.ax.contourf(self.nu3, self.nu1, self.dA, levels=levels, colors=colors)
            print(self.nu1)
            if self.lines==True:
                self.ax.contour(self.nu3, self.nu1, self.dA, colors='k',levels=levels[levels>0.06*self.scale[1]], linewidths=0.3, alpha=0.8, linestyles='-')
                self.ax.contour(self.nu3, self.nu1, self.dA, colors='k',levels=levels[levels<0.06*self.scale[0]], linewidths=0.3, alpha=0.8, linestyles='-')
        else:
            D = self.ax.contourf(self.nu1, self.nu3, self.dA.T, levels=levels, colors=colors)
            if self.lines==True:
                self.ax.contour(self.nu1, self.nu3, self.dA.T, colors='k',levels=levels[levels>0.06*self.scale[1]], linewidths=0.3, alpha=0.8, linestyles='-')
                self.ax.contour(self.nu1, self.nu3, self.dA.T, colors='k',levels=levels[levels<0.06*self.scale[0]], linewidths=0.3, alpha=0.8, linestyles='-')
        self.cbar = plt.colorbar(D, ax=self.ax)

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        if not self.norm==True:
            self.cbar.set_label(r'$\Delta{A} \, \times 10^{3}$')  
        else:
            self.cbar.set_label(r'norm. $\Delta{A}$')
        self.cbar.set_ticks([self.scale[0], self.scale[1]])

        self.cbar.set_ticks([self.scale[0], 0, self.scale[1]])
        self.cbar.set_label(r'$\Delta A$ / mOD$\cdot$cm$^{-1/2}$')
        if self.diagonal==True:
            self.ax.plot(self.nu1, self.nu1, color='k', lw=0.75, linestyle='-', alpha=0.6)
        if self.flip==False:
            self.ax.set_ylabel(r'$\tilde{\nu}_{3} = \tilde{\nu}_{\text{probe}} / \text{cm}^{-1}$')
            self.ax.set_xlabel(r'$\tilde{\nu}_{1} =\tilde{\nu}_{\text{pump}} / \text{cm}^{-1}$')
        else:
            self.ax.set_xlabel(r'$\tilde{\nu}_{3} = \tilde{\nu}_{\text{probe}} / \text{cm}^{-1}$')
            self.ax.set_ylabel(r'$\tilde{\nu}_{1} = \tilde{\nu}_{\text{pump}} / \text{cm}^{-1}$')            
        self.ax.set_title(r'$t_2 = %.2g$ ps'%(self.delay))

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
            self.fig.savefig('%s_contour.pdf'%(self.file[:self.file.find('.')]), transparent=True)
        plt.show()