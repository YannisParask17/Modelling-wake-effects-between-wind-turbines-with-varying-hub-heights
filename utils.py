import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# Latex font
# mpl.style.use('classic')
plt.rcParams['font.family']             = 'STIXGeneral'
plt.rcParams["legend.scatterpoints"]    = 2
plt.rcParams["legend.numpoints"]        = 1
plt.rcParams.update({'font.size': 12})
mpl.rcParams['lines.linewidth']         = 1.5


class wind_farm_layout:
    def __init__(self, nrows, ncols, Sx, Sy, Dref, zref) -> None:
        """
        Class that can generate box-shaped wind farm layouts with equidistant turbine spacing 
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

        
        
    def visualise_wf(self, normalised=True):
        """_summary_

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
        # plt.plot(wt_x_VS/D, wt_y_VS/D, 'or')
        plt.show()



if __name__ == '__main__':
    wf_layout = wind_farm_layout(nrows=4,ncols=4,Dref=90, zref=80, Sx=5, Sy=2.5)
    wf_layout.get_wt_locations(D_VS=172, z_VS=199)
    wf_layout.visualise_wf(normalised=False)
    print(wf_layout.type_i_pwe)
    print(wf_layout.wt_x_pwe)