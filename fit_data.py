
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

class Fit_data:
    """A class for data."""
    
    def __init__(self, x, y, xtitle, ytitle, name, xunits, yunits, log):
        self.xrange = [x.min(), x.max()]
        self.yrange = [y.min(), y.max()]
    
        self.x = x
        self.y = y
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.log = log
        self.name = name
        
    def func(self, x):
        if self.log[0] == 1:
            fit_x = np.log10(self.x)
        else: fit_x = self.x
        
        if self.log[1] == 1:
            fit_y = np.log10(self.y)
        else: fit_y = self.y
        
        func = interpolate.interp1d(fit_x, fit_y, kind = 'quadratic')    

        return func(x)

    def show(self):
        ax = plt.subplot(111)
    
        if self.log is not None:
            if self.log[0] == 1:
                ax.set_xscale('log')
            if self.log[1] == 1:
                ax.set_yscale('log')
    
        ax.set_ylabel(self.ytitle)
        ax.set_xlabel(self.xtitle)
        ax.set_title(self.name)
        
        ax.plot(self.x, self.y, "o")
        num_points = self.x.shape[0]

        fit_x = np.linspace(self.xrange[0], self.xrange[1], num = num_points*10)
        fit_y = self.func(fit_x)
        ax.plot(fit_x, fit_y, "-")
        
        plt.show()
        