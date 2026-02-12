import numpy as np
from matplotlib import pyplot as plt
from VautheyLab.miscellaneous import moving_average
from VautheyLab.miscellaneous import solvs
np.seterr(divide='ignore')
import os
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))

# Class to plot Emission spectra
class Emission:
    # initialize class
    def __init__(self, files_em, figsize=None, delimiter_em = ',',
                 cuts_em=[None], norm_em=True, ylim_em=None, ylim_abs=None, xlim=None, xticks=True,
                 yticks=True, fill=False, inv=True, TDM=False, colors_em = ['r'], files_abs = [None], norm_abs=True, baseline_abs=[None],
                 cuts_abs = [None], labels_abs=[''], normat_abs = [None], colors_abs = ['r'],
                 labels_em=[''], title='', units='wn', plot_now=True, corr=None,
                 export=False, MA_abs=[None], MA_npoints_abs=[10], MA_em=[None], MA_npoints_em=[10], tightlayout=True, savefig=False,
                 linestyle=None, baseline_em=[None], normat_em=[None], waterfall=None, waterfall_shift=1.5, outside=False, fontsize=10, legend=True,
                 files_ex = [None], cuts_ex = [None],  labels_ex=[''], normat_ex = [None], colors_ex = ['r'], norm_ex = True, 
                 MA_npoints_ex = [5], MA_ex = [None]):
        
        # set relevant parameters
        self.figsize = figsize

        # get emission files
        self.files_em = files_em
        self.corr = corr

        self.files_ex = files_ex

        # get wavelength cuts
        self.cuts_em = cuts_em
        if self.cuts_em[0]!=None:
            if len(self.cuts_em)==1 and len(self.files_em)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_em) - 1):
                    self.cuts_em.append(self.cuts_em[0])
        # decide whether to normalize or not
        self.norm_em = norm_em
        # to normalize at a specific wl
        self.normat_em = normat_em
        self.cuts_em = cuts_em
        if self.normat_em[0]==None:
            if len(self.normat_em)==1 and len(self.files_em)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_em) - 1):
                    self.normat_em.append(self.normat_em[0])
        self.baseline_em = baseline_em
        self.labels_em = labels_em

        # excitaion spectra
        self.colors_ex = colors_ex
        self.cuts_ex = cuts_ex
        if self.cuts_ex[0]!=None:
            if len(self.cuts_ex)==1 and len(self.files_ex)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_ex) - 1):
                    self.cuts_ex.append(self.cuts_ex[0])
        # decide whether to normalize or not
        self.norm_ex = norm_ex
        # to normalize at a specific wl
        self.normat_ex = normat_ex
        self.cuts_ex = cuts_ex
        if self.normat_ex[0]==None:
            if len(self.normat_ex)==1 and len(self.files_ex)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_ex) - 1):
                    self.normat_ex.append(self.normat_ex[0])
        self.labels_ex = labels_ex

        # get absorption files
        self.files_abs = files_abs
        # get wavelength cuts
        self.cuts_abs = cuts_abs
        if self.cuts_abs[0]!=None:
            if len(self.cuts_abs)==1 and len(self.files_abs)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_abs) - 1):
                    self.cuts_abs.append(self.cuts_abs[0])
        # decide whether to normalize or not
        self.norm_abs = norm_abs
        # to normalize at a specific wl
        self.normat_abs = normat_abs
        self.cuts_abs = cuts_abs
        if self.normat_abs[0]==None:
            if len(self.normat_abs)==1 and len(self.files_abs)>1:
                # use same range if only one cuts specified
                for i in range(len(self.files_abs) - 1):
                    self.normat_abs.append(self.normat_abs[0])
        self.baseline_abs = baseline_abs
        self.labels_abs = labels_abs

        # delimiter for files
        self.delimiter_em = delimiter_em

        # specify x and y limits for plot
        self.ylim_em = ylim_em
        self.ylim_abs = ylim_abs
        self.xlim = xlim
        # specify xticks
        self.xticks = xticks
        self.yticks = yticks
        # decide whether to have spectra filled
        self.fill = fill
        # specify whether you want to invert the x-axis when using kK or not
        self.inv = inv
        # transition dipole moment representation
        self.TDM = TDM
        # specify list of colors
        self.colors_abs = colors_abs
        self.colors_em = colors_em
        # sepcify title of your plot
        self.title = title
        # sepecify your units
        self.units = units
        self.outside = outside
        self.legend = legend
        self.fontsize = fontsize

        # whether to export data or not
        self.export = export

        # moving average?
        self.MA_abs = MA_abs
        self.MA_npoints_abs = MA_npoints_abs
        if self.MA_abs[0]!=None:
            if len(self.MA_abs)==1 and len(self.MA_abs)>1:
                # use same range if only one cuts specified
                for i in range(len(self.MA_abs) - 1):
                    self.MA_abs.append(self.MA_abs[0])
                    self.MA_npoints_abs.append(self.MA_npoints_abs[0])

        # moving average?
        self.MA_ex = MA_ex
        self.MA_npoints_ex = MA_npoints_ex
        if self.MA_ex[0]!=None:
            if len(self.MA_ex)==1 and len(self.MA_ex)>1:
                # use same range if only one cuts specified
                for i in range(len(self.MA_ex) - 1):
                    self.MA_ex.append(self.MA_ex[0])
                    self.MA_npoints_em.append(self.MA_npoints_ex[0])

        # moving average?
        self.MA_em = MA_em
        self.MA_npoints_em = MA_npoints_em
        if self.MA_em[0]!=None:
            if len(self.MA_em)==1 and len(self.MA_em)>1:
                # use same range if only one cuts specified
                for i in range(len(self.MA_em) - 1):
                    self.MA_em.append(self.MA_em[0])
                    self.MA_npoints_em.append(self.MA_npoints_em[0])

        # tight layout
        self.tightlayout = tightlayout
        self.savefig = savefig

        # for waterfall plot
        self.waterfall = waterfall
        self.waterfall_shift = waterfall_shift

        # line style 
        if linestyle==None:
            self.linestyle = ['-' for i in range(len(self.files_em))]
        else:
            self.linestyle = linestyle
        print(linestyle)

        # if no labels specified
        if self.labels_em==[''] and len(self.files_em)>1:
            self.labels_em = ['' for i in range(len(self.files_em))]

        # if only one baseline use same for all files
        if self.baseline_em!=None:
            if len(self.baseline_em)==1 and len(self.files_em)>1:
                self.baseline_em = [baseline_em[0] for i in range(len(self.files_em))]

        # if no labels specified
        if self.labels_abs==[''] and len(self.files_abs)>1:
            self.labels_abs = ['' for i in range(len(self.files_abs))]

        # if only one baseline use same for all files
        if self.baseline_abs!=None:
            if len(self.baseline_abs)==1 and len(self.files_abs)>1:
                self.baseline_abs = [baseline_abs[0] for i in range(len(self.files_abs))]

        # intiialize figure
        self.fig, self.ax = plt.subplots(1,1)
        # read data
        self.read_data_em()
        if self.files_abs[0]!=None:
            self.read_data_abs()
        if self.files_ex[0]!=None:
            self.read_data_ex()
        # plot data
        if plot_now==True:
            self.plot()

    # read in data
    def read_data_abs(self):
        # initialize list of absorbances
        self.A = []
        self.wl_abs = []
        self.wn_abs = []
        # go through all files
        for i in range(len(self.files_abs)):
            # read file
            if not '.txt' in self.files_abs[i]:
                data = np.loadtxt('%s'%self.files_abs[i], delimiter=',', usecols=[0,1], skiprows=2)
            else:
                data = np.loadtxt('%s'%self.files_abs[i], delimiter=',', usecols=[0,2], skiprows=1)
            
            # get wavelength
            if data[:,0][0]<100:
                self.wl_abs.append((1/data[:,0])*10**4)
                self.wn_abs.append((1/self.wl_abs[i])*10**4)
            else:
                self.wl_abs.append(data[:,0])
                self.wn_abs.append((1/data[:,0])*10**4)
            # get absorbance
            self.A.append(data[:,1]) 
            
            # correct baseline 
            if self.baseline_abs[0]!=None:
                if self.baseline_abs[i]!=None:
                    bdata = np.loadtxt(self.baseline_abs[i], delimiter=',', usecols=[0,1], skiprows=2)
                    self.A[i] = self.A[i] - bdata[:,1]

            # cut data 
            if self.cuts_abs[0]!=None:
                self.A[i] = self.A[i][(self.wl_abs[i]>self.cuts_abs[i][0])&(self.wl_abs[i]<self.cuts_abs[i][1])]
                self.wl_abs[i] = self.wl_abs[i][(self.wl_abs[i]>self.cuts_abs[i][0])&(self.wl_abs[i]<self.cuts_abs[i][1])]
                self.wn_abs[i] = (1/self.wl_abs[i])*10**4

            # normalize
            if self.norm_abs == True: 
                if self.normat_abs[i]==None:
                    self.A[i] = self.A[i]/np.max(self.A[i])
                else:
                    self.A[i] = self.A[i]/self.A[i][np.argmin(np.abs(self.wl_abs[i] - self.normat_abs[i]))]
            # export if wanted
            if self.export==True:
                if self.norm_abs == True: 
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  norm. absorbance'
                else:
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance'
                np.savetxt('abs_%s.txt'%(self.files_abs[i][:self.files_abs[i].find('.')]), 
                            np.column_stack([self.wl_abs[i], self.wn_abs[i], self.A[i]]),
                            header=head, delimiter=',')

    # read in data
    def read_data_em(self):
        # initialize list of absorbances
        self.I_wl = []
        self.I_wn = []
        self.wl_em = []
        self.wn_em = []
        # go through all files
        for i in range(len(self.files_em)):
            # read file
            if not '.txt' in self.files_em[i]:
                data = np.loadtxt('%s'%self.files_em[i], delimiter=self.delimiter_em, usecols=[0,1], skiprows=2)
            else:
                data = np.loadtxt('%s'%self.files_em[i], delimiter=',', usecols=[0,2], skiprows=1)
            
            # get wavelength
            if data[:,0][0]<100:
                self.wl_em.append((1/data[:,0])*10**4)
                self.wn_em.append((1/self.wl_em[i])*10**4)
            else:
                self.wl_em.append(data[:,0])
                self.wn_em.append((1/data[:,0])*10**4)

            # get emission
            self.I_wl.append(data[:,1]) 
            
            # correct baseline 
            if self.baseline_em[0]!=None:
                bdata = np.loadtxt(self.baseline_em[i], delimiter=',', usecols=[0,1], skiprows=2)
                self.I_wl[i] = self.I_wl[i] - bdata[:,1]

            # cut data 
            if self.cuts_em[0]!=None:
                self.I_wl[i] = self.I_wl[i][(self.wl_em[i]>self.cuts_em[i][0])&(self.wl_em[i]<self.cuts_em[i][1])]
                self.wl_em[i] = self.wl_em[i][(self.wl_em[i]>self.cuts_em[i][0])&(self.wl_em[i]<self.cuts_em[i][1])]
                self.wn_em[i] = (1/self.wl_em[i])*10**4
            
            # correct data 
            if self.corr==True:
                # subtract dark counts
                self.I_wl[i] = self.I_wl[i]
                corr = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'FMax_lamp_20151217.txt'))
                c = np.interp(self.wl_em[i], corr[:,0], corr[:,1])
                self.I_wl[i] = c*self.I_wl[i]

            # convert intensity to wavenumber
            self.I_wn.append(self.I_wl[i] * self.wl_em[i]**2)

            # normalize
            if self.norm_em == True: 
                if self.normat_em[i]==None:
                    self.I_wl[i] = self.I_wl[i]/np.max(self.I_wl[i])
                    self.I_wn[i] = self.I_wn[i]/np.max(self.I_wn[i])
                else:
                    self.I_wl[i] = self.I_wl[i]/self.I_wl[i][np.argmin(np.abs(self.wl_em[i] - self.normat_em[i]))]
                    self.I_wn[i] = self.I_wn[i]/self.I_wn[i][np.argmin(np.abs(self.wl_em[i] - self.normat_em[i]))]

            # export if wanted
            if self.export==True:
                if self.norm_em == True: 
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  norm. Intensity for wavelength, norm. Intensity for wavenumber'
                else:
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity for wavelength / a.u., Intensity for wavenumber / a.u.'
                if self.MA_em[i]==False:
                    np.savetxt('em_%s.txt'%(self.files_em[i][:self.files_em[i].find('.')]), 
                                np.column_stack([self.wl_em[i], self.wn_em[i], self.I_wl[i], self.I_wn[i]]),
                                header=head, delimiter=',')
                else:
                    np.savetxt('em_%s.txt'%(self.files_em[i][:self.files_em[i].find('.')]), 
                                np.column_stack([moving_average(self.wl_em[i], self.MA_npoints_em[i]), moving_average(self.wn_em[i], self.MA_npoints_em[i]), 
                                                 moving_average(self.I_wl[i], self.MA_npoints_em[i]), moving_average(self.I_wn[i], self.MA_npoints_em[i])]),
                                header=head, delimiter=',')                    
    
    # read in data
    def read_data_ex(self):
        # initialize list of absorbances
        self.wlex = []
        self.Iex = []
        self.wnex = []
        # go through all files
        for i in range(len(self.files_ex)):
            # read file
            data = np.loadtxt('%s'%self.files_ex[i], delimiter=',', usecols=[0,1], skiprows=2)
            self.wlex.append(data[:,0])
            self.wnex.append((1/data[:,0])*10**4)

            # get emission
            self.Iex.append(data[:,1]) 
            
            # cut data 
            if self.cuts_ex[0]!=None:
                self.Iex[i] = self.Iex[i][(self.wlex[i]>self.cuts_ex[i][0])&(self.wlex[i]<self.cuts_ex[i][1])]
                self.wlex[i] = self.wlex[i][(self.wlex[i]>self.cuts_ex[i][0])&(self.wlex[i]<self.cuts_ex[i][1])]
                self.wnex[i] = (1/self.wlex[i])*10**4
            
            # normalize
            if self.norm_ex == True: 
                if self.normat_ex[i]==None:
                    self.Iex[i] = self.Iex[i]/np.max(self.Iex[i])
                else:
                    self.Iex[i] = self.Iex[i]/self.Iex[i][np.argmin(np.abs(self.wlex[i] - self.normat_ex[i]))]

            # export if wanted
            if self.export==True:
                if self.norm_ex == True: 
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  norm. Intensity'
                else:
                    head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  Intensity / a.u.'
                np.savetxt('ex_%s.txt'%(self.files_ex[i][:self.files_ex[i].find('.')]), 
                            np.column_stack([self.wlex[i], self.wnex[i], self.Iex[i]]),
                            header=head, delimiter=',')
                
    # find maximum
    def find_max(self, file_index, which):
        if which=='abs':
            Amax = np.max(self.A[file_index])
            wlmax = self.wl_abs[file_index][self.A[file_index] == np.max(self.A[file_index])]
            wnmax = self.wn_abs[file_index][self.A[file_index] == np.max(self.A[file_index])]
            print("%s --> Amax = %.3g at wl = %.3g nm / wn = %.3g kK"%(self.files_abs[file_index], Amax, wlmax, wnmax))
            return Amax, wlmax, wnmax
        if which=='em':
            Imax_wl = np.max(self.I_wl[file_index])
            Imax_wn = np.max(self.I_wn[file_index])
            wlmax = self.wl_em[file_index][self.I_wl[file_index] == np.max(self.I_wl[file_index])]
            wnmax = self.wn_em[file_index][self.I_wn[file_index] == np.max(self.I_wn[file_index])]
            print("%s --> Imax = %.3g / %.3g at wl = %.3g nm / wn = %.3g kK"%(self.files_em[file_index], Imax_wl, Imax_wn, wlmax, wnmax))
            return Imax_wl, Imax_wn, wlmax, wnmax
    
    def plot(self):
        # for waterfall plot
        shift_abs = 0
        shift_em = 0
        # plot absorption spectra
        if self.files_abs[0]!=None:
            for i in range(len(self.A)):
                if self.units=='wl':
                    x = self.wl_abs[i]
                else:
                    x = self.wn_abs[i]

                # transition dipole moment representation
                if self.TDM==True:
                    y = self.A[i]/x
                    y = y/np.max(y)
                else:
                    y = self.A[i]

                # shift if waterfall
                if self.waterfall!=None:
                    y += shift_abs
                    shift_abs+=self.waterfall_shift

                # plot
                if self.MA_abs[0]==None:
                    self.ax.plot(x, y, color=self.colors_abs[i], label=self.labels_abs[i])
                    if self.fill == True:
                        self.ax.fill_between(x, np.min(y), y, color=self.colors_abs[i], alpha=0.1)  
                else:
                    self.ax.plot(x, y, '-', color=self.colors_abs[i], alpha=0.2)
                    self.ax.plot(moving_average(x, self.MA_npoints_abs), moving_average(y, self.MA_npoints_abs), '-', color=self.colors_abs[i], label=self.labels_abs[i])
                    if self.fill == True:
                        self.ax.fill_between(moving_average(x, self.MA_npoints_abs), np.min(y), moving_average(y, self.MA_npoints_abs), color=self.colors_abs[i], alpha=0.1)                  

        # plot excitation spectra
        if self.files_ex[0]!=None:
            for i in range(len(self.Iex)):
                if self.units=='wl':
                    x = self.wlex[i]
                else:
                    x = self.wnex[i]
                y = self.Iex[i]
                # plot
                if self.MA_ex[0]==None:
                    self.ax.plot(x, y, color=self.colors_ex[i], label=self.labels_ex[i], linestyle=self.linestyle[i])
                    if self.fill == True:
                        self.ax.fill_between(x, np.min(y), y, color=self.colors_ex[i], alpha=0.1)  
                else:
                    self.ax.plot(x, y, '-', color=self.colors_ex[i], alpha=0.2)
                    self.ax.plot(moving_average(x, self.MA_npoints_ex[i]), moving_average(y, self.MA_npoints_ex[i]), '-', color=self.colors_ex[i], label=self.labels_ex[i])
                    if self.fill == True:
                        self.ax.fill_between(moving_average(x, self.MA_npoints_ex[i]), np.min(y), moving_average(y, self.MA_npoints_ex[i]), color=self.colors_ex[i], alpha=0.1)                  

        # plot emission spectra
        for i in range(len(self.I_wl)):

            if self.units=='wl':
                x = self.wl_em[i]
                y = self.I_wl[i]
            else:
                x = self.wn_em[i]
                y = self.I_wn[i]

            # transition dipole moment representation
            if self.TDM==True:
                y = self.I_wn[i]/(self.wn_em[i]**3)
                y = y/np.max(y)

            # shift if waterfall
            if self.waterfall!=None:
                y += shift_em
                shift_em+=self.waterfall_shift

            # plot
            if self.MA_em[0]==None:
                self.ax.plot(x, y, color=self.colors_em[i], label=self.labels_em[i], linestyle=self.linestyle[i])
                if self.fill == True:
                    self.ax.fill_between(x, np.min(y), y, color=self.colors_em[i], alpha=0.1)  
            else:
                if self.norm_em==True:
                    norm = np.max(moving_average(y, self.MA_npoints_em[i]))
                else:
                    norm = 1
                self.ax.plot(x, y/norm, '-', color=self.colors_em[i], alpha=0.2)
                self.ax.plot(moving_average(x, self.MA_npoints_em[i]), moving_average(y, self.MA_npoints_em[i])/norm, '-', color=self.colors_em[i], label=self.labels_em[i])
                if self.fill == True:
                    self.ax.fill_between(moving_average(x, self.MA_npoints_em[i]), np.min(y), moving_average(y, self.MA_npoints_em[i])/norm, color=self.colors_em[i], alpha=0.1)  

    # function to plot absorption solvatochromism
    def solvchrom_abs(self, solv, lim=None, save=False):
        fig1, ax1 = plt.subplots(1,1,figsize=(5, 3.5))
        fig2, ax2 = plt.subplots(1,1,figsize=(5, 3.5))
        fig3, ax3 = plt.subplots(1,1,figsize=(5, 3.5))
        abs_max = []
        for i in range(len(solv)):
            n = solvs[solv[i]][0]
            er = solvs[solv[i]][1]
            fe = (2*(er - 1)/(2*er+1))
            fn = (2*(n**2 - 1)/(2*n**2+1))
            df = fe - fn
            abs_max.append(self.wn_abs[i][self.A[i]==np.max(self.A[i])])
            ax1.plot(df, self.wn_abs[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors_abs[i])
            ax2.plot(fe, self.wn_abs[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors_abs[i])
            ax3.plot(fn, self.wn_abs[i][self.A[i]==np.max(self.A[i])], 'o', color=self.colors_abs[i])
        np.savetxt("export_abs_solvchrom.txt",  abs_max, delimiter=",", header='wavenumber at maxiumum absorbance')
        ax1.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax2.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax3.set_ylabel(r'$\Tilde{\nu}_{\text{max, abs}} / 10^3$ cm$^{-1}$')
        ax1.set_xlabel(r'$\Delta f = f(\varepsilon_r) - f(n^2)$')
        ax2.set_xlabel(r'$f(\varepsilon_r)$')
        ax3.set_xlabel(r'$f(n^2)$')
        if lim!=None:
            ax1.set_ylim(lim)
            ax2.set_ylim(lim)
            ax3.set_ylim(lim)
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        if save!=False:
            fig1.savefig('df_abs.svg', transparent=True)
            fig2.savefig('f_epsr_abs.svg', transparent=True)
            fig3.savefig('f_n2_abs.svg', transparent=True)

    # function to plot emission solvatochromism
    def solvchrom_em(self, solv, lim=None, save=False):
        fig1, ax1 = plt.subplots(1,1,figsize=(5, 3.5))
        fig2, ax2 = plt.subplots(1,1,figsize=(5, 3.5))
        fig3, ax3 = plt.subplots(1,1,figsize=(5, 3.5))
        max_em = []  
        for i in range(len(solv)):
            n = solvs[solv[i]][0]
            er = solvs[solv[i]][1]
            fe = (2*(er - 1)/(2*er+1))
            fn = (2*(n**2 - 1)/(2*n**2+1))
            df = fe - fn
            max_em.append(self.wn_em[i][self.I_wn[i]==np.max(self.I_wn[i])])
            ax1.plot(df, self.wn_em[i][self.I_wn[i]==np.max(self.I_wn[i])], 'x', color=self.colors_em[i])
            ax2.plot(fe, self.wn_em[i][self.I_wn[i]==np.max(self.I_wn[i])], 'x', color=self.colors_em[i])
            ax3.plot(fn, self.wn_em[i][self.I_wn[i]==np.max(self.I_wn[i])], 'x', color=self.colors_em[i])
        np.savetxt("export_em_solvchrom.txt",  max_em, delimiter=",", header='wavenumber at maxiumum emission')
        ax1.set_ylabel(r'$\Tilde{\nu}_{\text{max, em}} / 10^3$ cm$^{-1}$')
        ax2.set_ylabel(r'$\Tilde{\nu}_{\text{max, em}} / 10^3$ cm$^{-1}$')
        ax3.set_ylabel(r'$\Tilde{\nu}_{\text{max, em}} / 10^3$ cm$^{-1}$')
        ax1.set_xlabel(r'$\Delta f = f(\varepsilon_r) - f(n^2)$')
        ax2.set_xlabel(r'$f(\varepsilon_r)$')
        ax3.set_xlabel(r'$f(n^2)$')
        if lim!=None:
            ax1.set_ylim(lim)
            ax2.set_ylim(lim)
            ax3.set_ylim(lim)
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        if save!=False:
            fig1.savefig('df_em.svg', transparent=True)
            fig2.savefig('f_epsr_em.svg', transparent=True)
            fig3.savefig('f_n2_em.svg', transparent=True)

    def show(self):
        if self.figsize!=None:
            self.fig.set_size_inches(self.figsize)
        
        if self.units=='wn':
            self.ax.set_xlabel(r'$\tilde{\nu} / 10^{3}\,\text{cm}^{-1}$')
            if self.inv==True:
                self.ax.invert_xaxis()    
            axsec = self.ax.secondary_xaxis('top', functions=(lambda x: (1/x)*10**4, lambda x: (1/x)*10**4))
            axsec.set_xlabel(r'$\lambda / $ nm') 

        if self.units=='wl':
            self.ax.set_xlabel(r'$\lambda / $ nm')

        if self.legend==True:
            if self.outside==True:
                self.ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fontsize)
            else:
                self.ax.legend(fontsize=self.fontsize)

        self.ax.axhline(y=0, color='k')

        if self.files_abs[0]!=None:
            self.ax2 = self.ax.twinx()

        if self.norm_em == True:
            if not self.TDM == True:
                if self.files_abs[0]!=None:
                    self.ax.set_ylabel('norm. Absorbance')
                    self.ax2.set_ylabel('norm. Intensity')
                else:
                    self.ax.set_ylabel('norm. Intensity')
            else:
                self.ax.set_ylabel(r'$\varepsilon(\tilde{\nu})/\tilde{\nu}$') 
                if self.files_abs[0]!=None:
                    self.ax2.set_ylabel(r'$F(\tilde{\nu})/\tilde{\nu}^3$') 
        else:
            if self.files_abs[0]==None:
                self.ax.set_ylabel('Intensity / a.u.')   
            else:
                self.ax.set_ylabel('Absorbance')
                self.ax2.set_ylabel('Intensity / a.u.')                         


        if self.yticks == False:
            self.ax.set_yticks([])
            if self.files_abs[0]!=None:
                self.ax2.set_yticks([])

        if self.ylim_em != None:
            if self.files_abs[0]==None:
                self.ax.set_ylim(self.ylim_em)
            else:
                self.ax.set_ylim(self.ylim_abs)
                self.ax2.set_ylim(self.ylim_em)

        if self.xlim != None: 
            self.ax.set_xlim(self.xlim)
        if self.xticks == False:
            self.ax.set_xticks([])

        if self.title != None:
            self.ax.set_title(self.title)
        if self.tightlayout == True:
            self.fig.tight_layout()
        if self.savefig == True: 
            fname = input("Filename? : ")
            self.fig.savefig('%s.svg'%(fname), transparent=True)
        plt.show()