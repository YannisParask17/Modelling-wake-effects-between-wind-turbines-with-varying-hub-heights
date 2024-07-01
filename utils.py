import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
import xarray


# Latex font
mpl.style.use('classic')
plt.rcParams['font.family']             = 'STIXGeneral'
plt.rcParams["legend.scatterpoints"]    = 2
plt.rcParams["legend.numpoints"]        = 1
plt.rcParams.update({'font.size': 12})
mpl.rcParams['lines.linewidth']         = 1.5


class wind_farm_layout:
    def __init__(self, nrows, ncols, Sx, Sy, Dref, zref) -> None:
        """
        Class that can generate square-shaped wind farm layouts with equidistant turbine spacing 
        in the streamwise and spanwise directions defined by the user. 
        Option for vertically-staggered arrangements.

        Parameters
        ----------
        nrows : int
            number of turbines in the streamwise direction
        ncols : int
            number of turbines in the spanwise direction
        Sx : float
            Streamwise turbine interspacing in rotor diameters [Dref]
        Sy : float
            Spanwise turbine interspacing in rotor diameters [Dref]
        Dref : float
            reference rotor diameter. In vertically-staggered setups Dref corresponds to the short (baseline) turbines [m]
        zref : float
            reference hub height. In vertically-staggered setups zref corresponds to the short (baseline) turbines [m]
        """
        self.nrows  = nrows
        self.ncols  = ncols
        self.Sx     = Sx
        self.Sy     = Sy
        self.Dref   = Dref
        self.zref   = zref

    def get_wt_locations(self, D_VS=False, z_VS=False):
        """
        Output wind turbine sites with option for vertically-staggered setup. 
        The x (streamwise) and y(spanwise) locations and the turbine type indicators are stored
        in a form compatible to PyWakeEllipSys.  
        Note that a box layout is implemented for simplicity with constant spanwise and streamwise spacing  

        Parameters
        ----------
        D_VS : float, optional
            rotor diameter of the tall turbines in a vertically-staggered setup, by default False
        z_VS : float, optional
            hub height of the tall turbines in a vertically-staggered setup, by default False
        """

        wt_x_norm = np.arange(0, self.nrows*self.Sx, self.Sx)
        wt_y_norm = np.arange(0, self.ncols*self.Sy, self.Sy)

        wt_x        =   []
        wt_y        =   []
        # append x and y location of each wt (normalised by rotor diameter)
        for y in wt_y_norm:
            for x in wt_x_norm:
                wt_x.append(x)
                wt_y.append(y)
        type_i = np.zeros(len(wt_x))        
        wt_x = np.array(wt_x)
        wt_y = np.array(wt_y)

        # add tall turbine locations for vertically-staggered setup
        wt_x_VS     =   []
        wt_y_VS     =   [] 
        if D_VS and z_VS:
            for y in wt_y_norm:
                for x in wt_x_norm[:-1]:
                    wt_x_VS.append(x+self.Sx/2)
                    wt_y_VS.append(y)
        wt_x_VS     = np.array(wt_x_VS)
        wt_y_VS     = np.array(wt_y_VS)
        type_i_VS   = np.ones(len(wt_x_VS)) 

        
        # save variables in form compatible to PyWakeEllipSys
        self.wt_x       =   wt_x * self.Dref # rescale with rotor diameter (of short turbines)
        self.wt_y       =   wt_y * self.Dref # rescale with rotor diameter (of short turbines)
        if D_VS and z_VS:
            self.D_VS       =   D_VS
            self.z_VS       =   z_VS
            self.wt_x_VS    =   wt_x_VS * self.Dref   
            self.wt_y_VS    =   wt_y_VS * self.Dref
            self.wt_x_pwe   =   np.concatenate((self.wt_x, self.wt_x_VS))
            self.wt_y_pwe   =   np.concatenate((self.wt_y, self.wt_y_VS))
            self.type_i_pwe =   np.concatenate((type_i, type_i_VS))
        else:
            self.wt_x_pwe   =   self.wt_x
            self.wt_y_pwe   =   self.wt_y
            self.type_i_pwe =   type_i
    

    def center_wf(self):
        """
        Method that computes the turbine location bringing the origin (0,0) to
        the geometrical centre of the wind farm  
        """
        self.wt_x_pwe = self.wt_x_pwe - (self.wt_x_pwe.min() + self.wt_x_pwe.max())/2
        self.wt_y_pwe = self.wt_y_pwe - (self.wt_y_pwe.min() + self.wt_y_pwe.max())/2

        
        
    def visualise_wf(self, normalised=True, savename=False):
        """
        Visualise the generated wind farm layout 

        Parameters
        ----------
        normalised : bool, optional
            Normalise arrangement with Dref, by default True
        """

        if normalised:
            normalise   = self.Dref
            xlabel      = r'Streamwise coordinate $x/D$ [-]'
            ylabel      = r'Spanwise coordinate $y/D$ [-]'
        else:
            normalise   = 1
            xlabel      = r'Streamwise coordinate $x$ [m]'
            ylabel      = r'Spanwise coordinate $y$ [m]'

        plt.figure
        plt.plot(self.wt_x/normalise, self.wt_y/normalise, 'ob')
        if hasattr(self, 'D_VS') and hasattr(self, 'z_VS'): 
            plt.plot(self.wt_x_VS/normalise, self.wt_y_VS/normalise, 'or')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if savename:
            plt.savefig(savename, bbox_inches='tight')
        plt.show()



class precursor_post:
    def __init__(self) -> None:
        """
        Class to post process 1D ABL profiles produced 
        through precursor simulations in EllipSys1D
        """
        self.prec_results = {}

        
    def load_precursor(self, infile, case, Uinf=1.0):
        """
        Reads precursor result file and stores variables 

        Parameters
        ----------
        infile : string
            precursor result file location 
        case : string
            name of the case to be loaded
        Uinf : float, optional
            wind speed at hub height [m/s], by default 1.0 corresponding to a non-dimensional profile 
        """

        data = np.loadtxt(infile)
        z,u,v,nuT,tke,_,theta = data[:,0],  data[:,1],  data[:,2], data[:,3], data[:,4],  data[:,5], data[:,6]
        
        Umag    =   np.sqrt(u**2 + v**2)            # compute velocity magnitude
        phi     =   np.rad2deg(np.arctan(v/u))      # compute wind direction
        TI      =   np.sqrt(2/3 * tke) / Uinf       # compute turbulence intensity
        theta_0 =   288.5

        # storing variables
        self.prec_results[case] = {'z'     :   z,
                             'U'     :   u,
                             'V'     :   v,
                             'Umag'  :   Umag,
                             'phi'   :   phi,
                             'nuT'   :   nuT,
                             'TI'    :   TI} 


    def plot_precursor(self, case_list, variables=['Umag', 'TI', 'phi', 'nuT'], Dref=False, zref=False, savename=False, **kwargs):
        """
        Visualises 1D ABL profiles produced 
        through precursor simulations in EllipSys1D

        Parameters
        ----------
        case_list : list of case names to be plotted
            names of cases to be plotted
        variables : list, optional
            variables to be plotted, by default ['Umag', 'TI', 'phi', 'nuT']
        Dref : float, optional
            reference rotor diameter [m], by default False
        zref : float, optional
            reference hub height [m], by default False 
        savename : string, optional
            name of the file (with file extension) for saving, by default False
        """

        fig, axs = plt.subplots(1,len(variables),sharey=True, figsize=(9,2.5))

        for case in case_list:
            for ivar, var in enumerate(variables):
                if var == 'Umag':
                    xlabel = r'$\sqrt{\overline{u}^2+\overline{v}^2} / U_{\infty}$ [-]'
                elif var == 'TI':
                    xlabel = r'$\sqrt{2/3 k}/U_{ref}$ [-]'
                elif var == 'phi':
                    xlabel = r'$\overline{\phi}$ [deg]'
                elif var == 'nuT':
                    xlabel = r'$\nu_T$ [m$^2$s$^{-1}$]'
                
                axs[ivar].plot(self.prec_results[case][var], self.prec_results[case]['z'])
                axs[ivar].set_xlabel(xlabel)

        axs[0].set_ylabel(r'$z$ [m]')
        axs[-1].set_ylim(0,2000)

        if Dref and zref:
            rotor_upper = zref + Dref/2
            rotor_lower = zref - Dref/2
            [ax.axhline(y=rotor_upper, color='k', linestyle='--') for ax in axs]
            [ax.axhline(y=rotor_lower, color='k', linestyle='--') for ax in axs]

        [ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='both')) for ax in axs]

        if kwargs.get('labels'):
            ncol = len(kwargs['labels'])
            fig.legend(kwargs['labels'], ncol=ncol, loc='upper center', bbox_to_anchor=(0.5, 1.15))

        if savename:
            plt.savefig(savename,bbox_inches='tight')

        plt.show()

    # theta = theta_0 + N_b[case_idx]**2 * theta_0 / 9.81 * z
    



class windfarm_post():
    def __init__(self) -> None:
        """
        Class to post process 3D wind farm flow fields produced 
        through numerical simulations in EllipSys3D
        """
        self.wf_results = {}
        pass


    def load_wf_results(self, wind_farm, case, forcemethod, turbmodel, cells_D, Ti, Uinf, Dref):
        """_summary_

        Parameters
        ----------
        wind_farm : object
            wind farm object (see above wind_farm_layout class)
        case : string
            casename
        forcemethod : string
            Actuator disk force method
        turbmodel : string
            turbulence model
        cells_D : float
            wake refinement region resolution
        Ti : float
            turbulence intensity [-]
        Uinf : float
            wind speed at hub height [m/s]
        Dref : float
            reference rotor diameter [m]
        """

        wt_x    = wind_farm.wt_x
        wt_y    = wind_farm.wt_y
        # wt_x_VS = wind_farm.wt_x_VS
        # wt_y_VS = wind_farm.wt_y_VS
        wt_x_VS = False
        wt_y_VS = False

        folder  = f'run_{case}_{forcemethod}_{turbmodel}_{cells_D}cD_Ti{Ti}_zeta0_wdGrid270/post_flow_wd270_ws10/'
        infile  = folder + 'flowdata.nc'
        data    = xarray.open_dataset(infile)
        flowvarnames = list(data.keys())

        self.wf_results[case] = {'Uinf'         : Uinf,
                                 'TI'           : Ti,
                                 'forcemethod'  : forcemethod,
                                 'turbmodel'    : turbmodel,
                                 'cells_D'      : cells_D,
                                 'wt_x'         : wt_x,
                                 'wt_y'         : wt_y,
                                 'wt_x_VS'      : wt_x_VS,
                                 'wt_y_VS'      : wt_y_VS,
                                 'Dref'         : Dref,
                                 'data'         : data}
        
        # for ivar, var in enumerate(flowvarnames):
        #     self.wf_results[case][var] = data.variables[var]
        
    def compute_spanwise_average(self, case):
        """
        Computing spanwise average flow field over the wind farm width

        Parameters
        ----------
        case : string
            casename
        """

        # 
        data            = self.wf_results[case]['data']
        flowvarnames    = list(data.keys())

        mask = (data.y >= self.wf_results[case]['wt_y'].min() - self.wf_results[case]['Dref']/2) & (data.y <= self.wf_results[case]['wt_y'].max() + self.wf_results[case]['Dref']/2)
        print('computing spanwise average flow field over wind farm width')

        if 'Umag' not in data.variables:
            data['Umag'] = np.sqrt(data.U**2 + data.V**2)

        for ivar, var in enumerate(flowvarnames):
            print('computing average: ', var)
            data[f'{var}_avg']    = np.mean(data[var].where(mask, drop=True), axis = 1)

        self.wf_results[case]['data']   = data 
        


    def examine_dampinglayer(self, case):
        """
        Visualising (spanwise averaged over wind farm width) eddy viscosity field 
        in the x-z plane to check the effectiveness of the damping layer combined 
        with the ABL-N model (see van der Laan et al. 2024 https://doi.org/10.5194/wes-2024-23)

        Parameters
        ----------
        case : string
            casename
        """
        
        cmap        = mpl.cm.get_cmap('jet')   
        color_turb  = 'black'                     
        xlabel      = r'$x/D$ [-]'
        ylabel      = r'$z$ [m]'
        dashes      = (5,4)
        xlim        = [-600, 600]
        # xlim        = False
        ylim        = [0,4000]
        cbar_bottom = 0.11
        cbar_width  = 0.02
        cbar_x      = 0.91



        data        =   self.wf_results[case]['data']
        nuT_avg     =   data.muT_avg / 1.225
        nuT_norm    =   nuT_avg/1.784060e-5      

        X,Y = np.meshgrid(data.x/self.wf_results[case]['Dref'], data.z, indexing='ij')
        fig, ax = plt.subplots(figsize=(9,2.5))
        im1 = ax.pcolormesh(X,Y,nuT_norm,cmap=cmap,
                                  shading='gouraud',vmin=nuT_norm.min(),vmax=nuT_norm.max(),
                                  rasterized=True)     

        xloc = self.wf_results[case]['wt_x']
        yloc = self.wf_results[case]['wt_y']
        Ntx  = len(xloc)
        Nty  = len(yloc)
        zh = 80
        Dref = self.wf_results[case]['Dref']
        
        for turb in range(int(Ntx)):
            x_turb = np.array([xloc[turb]-(xloc.min() + xloc.max())/2, xloc[turb]-(xloc.min() + xloc.max())/2]) / Dref
            z_turb = np.array([zh-Dref/2, zh+Dref/2])

            ax.plot(x_turb, z_turb,
                        color=color_turb,
                        linewidth=1.25)   
        ax.plot([-600,600],[2000, 2000], '--k')
        ax.plot([50 + wf_control.wt_x_pwe.max()/Dref,50 + wf_control.wt_x_pwe.max()/Dref],[0, 2000], '--k')
        ax.set(ylim=ylim, xlim=xlim)
        cbar    = fig.add_axes([cbar_x,cbar_bottom,cbar_width,1-cbar_bottom*2-0.018])
        cb      = fig.colorbar(im1, cax=cbar)
        cb.set_label(r'$\nu_T/\nu$ [-]',fontsize=15)
        cb.formatter.set_powerlimits((0, 0))
        cb.formatter.set_scientific(True)
        # cb1.locator  = tick_locator
        cb.update_ticks()

        fig.savefig(f'{case}_damping_layer.pdf' ,bbox_inches='tight')
        plt.show()

        
    




if __name__ == '__main__':

    test = {'layout':           True,
            'precursor' :       True,
            'damping_layer':    True}



    if test['layout']:
        """
        Test wind_farm_layout() class
        """
        wf_control  =   wind_farm_layout(nrows=4, ncols=4, Dref=90, zref=80, Sx=5, Sy=2.5)
        wf_control.get_wt_locations()
        wf_control.center_wf()
        wf_vs       =   wind_farm_layout(nrows=4, ncols=4, Dref=90, zref=80, Sx=5, Sy=2.5)
        wf_vs.get_wt_locations(D_VS=90, z_VS=90)
        wf_vs.center_wf()
        print(wf_control.wt_x_pwe)    
        print(wf_control.wt_y_pwe)    
        print(wf_vs.wt_x_pwe)
        print(wf_vs.wt_y_pwe)    

        show = False
        if show:
            fig, axs = plt.subplots(1,2,figsize=(9,2.5),sharex=True,sharey=True)
            xlim = [-1,16]
            ylim = [-1,22]
            normalise = wf_control.Dref
            axs[0].plot(wf_control.wt_x/normalise, wf_control.wt_y/normalise, 'ob', label='V90')
            axs[0].set(title='Control Wind Farm', 
                    xlabel=r'Streamwise coordinate $x/D$ [-]', 
                    ylabel=r'Spanwise coordinate $y/D$ [-]',
                    xlim=xlim,
                    ylim=ylim)

            axs[1].plot(wf_control.wt_x/normalise, wf_control.wt_y/normalise, 'ob', label='V90')
            axs[1].plot(wf_vs.wt_x_VS/normalise, wf_vs.wt_y_VS/normalise, 'or', label='V172')
            axs[1].set(title='VS Wind Farm', xlabel=r'Streamwise coordinate $x/D$ [-]')

            plt.legend(ncol=2)
            plt.savefig('RQ2_windfarm.pdf', bbox_inches='tight')


    if test['precursor']: 
        """
        Test precursor_post() class
        """       
        
        z0_list     = np.array([0.056, 0.08, 0.1])
        Gref_list   = np.array([1.347707558274214e+01, 1.415954496656181e+01, 1.472063638999129e+01])
        Ro0_list    = Gref_list / (1.14e-4 * z0_list)
        folder_list = ['run_precursor_opt_cori_Nb_Ti0.1_z0_0.056_Ro0_2.11107e+06/', 
                       'run_precursor_opt_cori_Nb_Ti0.1_z0_0.08_Ro0_1.55258e+06/',
                       'run_precursor_opt_cori_Nb_Ti0.1_z0_0.1_Ro0_1.29128e+06/']
        cases = [f'$z_0 = {z0}m$' for z0 in z0_list]

        precursor   = precursor_post()
        for i, folder in enumerate(folder_list):
            infile = folder + 'result_ws1.dat'
            precursor.load_precursor(infile, cases[i])

        precursor.plot_precursor(case_list=cases,Dref=90,zref=80,savename='precursor_VS.pdf', **{'labels' : cases})
            
        plt.show()
    
    if test['damping_layer']:
        """
        Test windfarm_post() class
        """       
        case        =   'control'
        forcemethod =   'Fix_con'
        turbmodel   =   'keABLcnfP'
        cells_D     =   4
        Ti          =   0.1
        Uinf        =   1
        wind_farm   =   wf_control
        Dref        =   80

        results = windfarm_post()
        results.load_wf_results(wind_farm, case, forcemethod, turbmodel, cells_D, Ti, Uinf, Dref)

        print(results.wf_results[case]['wt_x'])
        print(results.wf_results[case]['wt_y'])
        
        print(results.wf_results[case]['data'].x.min(), results.wf_results[case]['data'].x.max())
        print(results.wf_results[case]['data'].y.min(), results.wf_results[case]['data'].y.max())

        results.compute_spanwise_average(case)
        results.examine_dampinglayer(case)
        


