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
    def __init__(self, nrows, ncols, Dref, Sx, Sy) -> None:
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
        Dref : float
            reference rotor diameter. In vertically-staggered setups Dref corresponds to the short (baseline) turbines [m]
        Sx : int
            Streamwise turbine interspacing in rotor diameters [Dref]
        Sy : int
            Spanwise turbine interspacing in rotor diameters [Dref]
        """
        self.nrows  = nrows
        self.ncols  = ncols
        self.Dref   = Dref
        self.Sx     = Sx
        self.Sy     = Sy

    def get_wt_locations(self, vertically_staggered=False):

        wt_x_norm = np.arange(0, self.nrows*self.Sx, self.Sx)
        wt_y_norm = np.arange(0, self.ncols*self.Sy, self.Sy)

        wt_x = []
        wt_y = []
        # append x and y location of each wt (normalised by rotor diameter)
        for y in wt_y_norm:
            for x in wt_x_norm:
                wt_x.append(x)
                wt_y.append(y)
        
        # rescale with rotor diameter (of short turbines)
        wt_x = np.array(wt_x)
        wt_y = np.array(wt_y)
        self.wt_x   =   wt_x * self.Dref 
        self.wt_y   =   wt_y * self.Dref

        wt_x_VS = []
        wt_y_VS = [] 
        if vertically_staggered:
            for y in wt_y_norm:
                for x in wt_x_norm[:-1]:
                    wt_x_VS.append(x+self.Sx/2)
                    wt_y_VS.append(y)
        wt_x_VS = np.array(wt_x_VS)
        wt_y_VS = np.array(wt_y_VS)
        

        self.wt_x       =   wt_x * self.Dref 
        self.wt_y       =   wt_y * self.Dref
        self.wt_x_VS    =   wt_x_VS * self.Dref   
        self.wt_y_VS    =   wt_y_VS * self.Dref

        
        
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
        plt.plot(self.wt_x_VS/normalise, self.wt_y_VS/normalise, 'or')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        # plt.plot(wt_x_VS/D, wt_y_VS/D, 'or')
        plt.show()



if __name__ == '__main__':
    wf_layout = wind_farm_layout(4,4,90, Sx=5, Sy=2.5)
    wf_layout.get_wt_locations(True)
    wf_layout.visualise_wf(normalised=False)
