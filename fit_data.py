from __future__ import absolute_import

__all__ = ["Fit_data"]

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
        self.xunits = xunits
        self.yunits = yunits
        
    def func(self, x):
        if self.log[0] == 1:
            fit_x = np.log10(self.x)
        else: fit_x = self.x
        
        if self.log[1] == 1:
            fit_y = np.log10(self.y)
            fill_value = -100
        else: 
            fit_y = self.y
            fill_value = 0
        
        f = interpolate.interp1d(fit_x, fit_y, kind = 3, bounds_error=False, fill_value = fill_value)    
        x_in = x
        if self.log[0] == 1:
            x_in = 10 ** x_in
        
        if self.log[1] == 1:
            f1 = lambda y: 10 ** f(y)
        else:
            f1 = f
        
        return f1(x_in)

    def show(self):
        ax = plt.subplot(111)
    
        if self.log is not None:
            if self.log[0] == 1:
                ax.set_xscale('log')
            if self.log[1] == 1:
                ax.set_yscale('log')
    
        ax.set_ylabel(self.ytitle + ' [' + self.yunits + ']')
        ax.set_xlabel(self.xtitle + ' [' + self.xunits + ']')
        ax.set_title(self.name)
        
        num_points = self.x.shape[0]

        fit_x = np.linspace(self.xrange[0], self.xrange[1], num = num_points*10)
        fit_y = self.func(fit_x)
        ax.plot(fit_x, fit_y, "-", color = 'blue')
        
        ax.plot(self.x, self.y, "o", color = 'red')

        plt.show()
