import numpy as np
import matplotlib.pyplot as plt
import sys
import platform
import datetime
from scipy import constants as con
from scipy.special import kv
import urllib2 as url
from bs4 import BeautifulSoup
import os

from scipy.integrate import quad
from scipy import interpolate


#data_dir = os.path.join(os.path.dirname(heroes.__file__), "util", "data") 
data_dir = '/Users/schriste/Dropbox/python/heroes/util/data/'


# densities 
# source; wolframalpha
density = {"air stp": 0.001204, "si": 2.33, "be": 1.848, "water": 1, "cadmium telluride": 6.2, 
"cesium iodide": 4.51, "gallium arsenide": 5.31, "mercuric iodide": 6.36, "lead glass": 6.22}

'''The X-ray transmission data comes from NIST 
	(http://www.nist.gov/pml/data/xraycoef/index.cfm)'''    

def xray_transmission(path_length_m, energy_kev, material='air stp'):
	"""Provide the X-ray transmission (0 to 1) in given a path length in meters at 
	a particular energy given in keV through a material with a constant density."""
		
	coefficients = mass_attenuation_coefficicent(energy_kev, material=material)
	transmission = np.exp(-coefficients * density_cgs.get(material) * path_length_m * 100.0)
	
	return transmission

def load_mass_attenuation_coefficients(material='air_dry_near_sea_level'):
	'''Load the mass attenuation coefficients (cm2/g) and mass energy-absorption coefficients (cm2/g)
	from the data files as a function of energy (MeV). The allowed materials are listed in density.'''

	filename = data_dir + 'XrayMassCoef_' + material.replace(' ', '_').capitalize() + '.txt'
	data = np.genfromtxt(filename, comments = ';', missing_values = ' ', skip_header = 8)
	
	return data

def mass_attenuation_coefficicent(energy_kev, material):
    """Returns the mass attenuation coefficient at an energy given in keV"""

    data = load_mass_attenuation_coefficients(material)

    # data is better behaved in log space
    data_energy_kev = np.log10(data[:,0]*1000)
    data_attenuation_coeff = np.log10(data[:,1])
    f = interpolate.interp1d(data_energy_kev, data_attenuation_coeff)    

    return 10 ** f(np.log10(energy_kev))

def plot_mass_attenuation_coefficient(material='air_dry_near_sea_level'):
	'''Plot the mass the mass attenuation coefficients and mass energy-absorption 
	coefficients for a named material. See load_mass_attenuation_coefficients definition
	for list of allowed materials.'''
	
	data = load_mass_attenuation_coefficients(material=material)
		
	energy_kev = data[:,0]
	mass_atten_coeff = data[:,1]
	mass_energy_atten_coeff = data[:,2]
	
	ax = plt.subplot(111)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Energy [keV]')
	ax.set_title(material.replace('_', ' ').capitalize())
	ax.set_ylabel(r'Mass Attenuation Coefficient [cm$^2$/g]')
		
	ax.plot(energy_kev, mass_atten_coeff)
	ax.plot(energy_kev, mass_energy_atten_coeff)

	ax.legend((r'$\mu/\rho$', r'$\mu_{en}/\rho$'))

	plt.show()
	
def xray_absorption(energy_kev, thickness_um, material='si'):
	'''Calculate the xray absorption in a material with a thickess (given in microns).'''
	
	return 1-xray_transmission(energy_kev, thickness_um/1e6, material=material)

def detector_efficiency(energy_kev, thickness_um, material='si'):
	'''Calculate the detector quantum efficiency (in percent) at a given energy'''
	
	return xray_absorption(energy_kev, thickness_um, material=material)*100.0

def load_attenuation_length(material='si'):
	filename = data_dir + material + '_xray_atten_length.txt'
	data = np.genfromtxt(filename, comments = ';', missing_values = ' ', skip_header = 3)
	
	return data

def xyplot(x, y, ytitle = None, xtitle = None, title = None, log = None):
	
	ax = plt.subplot(111)

	if log is not None:
		if log[0] == 1:
			ax.set_xscale('log')
		if log[1] == 1:
			ax.set_yscale('log')

	if ytitle is not None:
		ax.set_ylabel(ytitle)
	if xtitle is not None:
		ax.set_xlabel(xtitle)
	if title is not None:		
		ax.set_title(title)
		
	ax.plot(x, y)
	plt.show()
	
    
def thermal_bremsstrahlung_thin(energy_kev, kt):
    """This function calculates the optically thin continuum thermal bremsstrahlung
    photon flux incident on the Earth from an isothermal plasma on the Sun.  
    Normalization is for an emission measure on the Sun of 1.e49 cm-3
    
    function brem_49,e,kt , verbose=verbose

if keyword_set(verbose) then print, 'Differential Bremsstrahlung '+$
  'spectrum at Earth for emission measure of 1.e49.'
;
kt0 =( kt(0) > 0.1) ;protect against vectors for kt
result = (1.e8/9.26) * float(acgaunt(12.3985/E, KT0/.08617)) *exp(-(E/KT0 < 50)) /E / KT0^.5

return, result(*)
    
    """
    
    # kt0 =( kt(0) > 0.1) ; protect against vectors for kt
    result = (1.e8/9.26) * gaunt_factor(energy_kev, kt) * 1/(energy_kev * np.sqrt(kt)) * n.exp(- (energy_kev / kt))

def gaunt_factor(energy_kev, kt):
    a = 0.5*energy_kev/kt
    return np.exp(a)*kv(a)

def effective_area(energy_kev):
    """Returns the HEROES effective area in cm^2 at a particular energy given in keV."""
    
    data_energy_kev = np.arange(20,80,10)
    data_effective_area = np.array([80,75,60,40,15,5])
    f = interpolate.interp1d(data_energy_kev, data_effective_area)

    return f(energy_kev)

def detector_background(energy_kev):
    
    data_energy_kev = np.arange(20,80,10)
    data_det_background = np.array([2,2,2.5,3,3,3]) * 0.001
    f = interpolate.interp1d(data_energy_kev, data_det_background)

    return f(energy_kev)

def atmo_transmission(energy_kev):
    
    data_energy_kev = np.arange(20,80,10)
    data_atmo_transmission = np.array([0.26, 2.0, 3.2, 3.7, 4.2, 4.5]) * 0.1
    f = interpolate.interp1d(data_energy_kev, data_atmo_transmission)

    return f(energy_kev)

def detector_efficiency(energy_kev):

    data_energy_kev = np.arange(20,80,10)
    data_detector_efficiency = np.array([9.8, 9.2, 9.9, 9.7, 8.9, 7.7]) * 0.1

    f = interpolate.interp1d(data_energy_kev, data_detector_efficiency)

    return f(energy_kev)

def sensitivity(background_counts, flux_to_counts_conversion, statistical_significance=5):
    """Calculates the sensitivity of an instrument using the following formula
    
        K = signal / sqrt(signal + background)
        
    where K is the significance (in sigma). This equation solves to 
    
        Sensitivity Flux limit = (K^2 + sqrt(K^4 - 4*background)) / 2 * flux_to_counts_conversion
    
    """
    
    result = 1/(2 * flux_to_counts_conversion) * statistical_significance ** 2 + np.sqrt( statistical_significance ** 4 - 4 * background_counts )
    
    return result
    

def sensitivity(integration_time, de = 5, statistical_sig = 5):
    """Returns the HEROES sensitivity at a particular energy given in keV. 
    de is the width of the energy interval in keV"""
    
    energy_kev = np.arange(20,80,10)
    det_eff = detector_background(energy_kev)
    det_background = detector_background(energy_kev)
    eff_area = effective_area(energy_kev)
    
    det_efficiency = detector_efficiency(energy_kev)
    transmission = atmo_transmission(energy_kev)
    background_area = 8 * 0.04
    fraction_flux = 0.8
    
    a = statistical_sig ** 2 + np.sqrt(statistical_sig ** 4 + 4*statistical_sig ** 2 *
        det_background * de * background_area * integration_time)
        
    b = 2 * eff_area * de * integration_time * transmission * det_eff * fraction_flux
    
    return  a/b

def get_msis_atmosphere_density(latitude=55, longitude=45, height=100, start=0, stop=1000, step=10):
    vars = [5,11] # 5 is height, 11 is density g/cm^3
    addr = 'http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi'
    data = u'model=msis&year=2000&month=01&day=01&time_flag=0&hour=1.5&geo_flag=0.&latitude=55.&longitude=45.'
    data = data + u'&height=100.&profile=1&start=0.&stop=1000.&step=50.&f10_7=&f10_7_3=&ap=&format=0&'
    data = data + 'vars=0' + str(vars[0]) + '&vars=0' + str(vars[1])
    a = url.Request(addr, data)
    f = url.urlopen(a)

    f = url.urlopen(a)
    
    data = np.genfromtxt(f, skip_header = 18, skip_footer = 16)

    return data

def atmosphere_density(height_km, date = '2012/09/01 00:00', latitude=55, longitude=45, data = None):
    '''
    Returns the atmospheric density (in g/cm^-3) at a specific height (given in km)
    
    Source
    ------
    http://omniweb.gsfc.nasa.gov/vitmo/msis_vitmo.html
    '''
    
    # get the data from online
    if data is None:
        data = get_msis_atmosphere_density(date=date, latitude=latitude, longitude=longitude)
    
    f = interpolate.interp1d(data[:,0], data[:,1])
        
    return f(height_km)
    
def xray_transmission_in_atmosphere(height_km, view_angle=90, data = None):
    """Find the total mass of atmosphere above a height given in km"""
    
    # get the data from online
    if data is None:
        data = get_msis_atmosphere_density()
 
    data_height_cm = data[:,0]*100000
    data_density_cgs = data[:,1]/np.sin(np.deg2rad(view_angle))
    
    co = interpolate.interp1d(height_cm, data[:,1])
        
    coefficients = mass_attenuation_coefficicent(energy_kev, material=material)
    transmission = np.exp(-coefficients * density_cgs.get(material) * path_length_m * 100.0)

    return quad(f, height_km, data[:,0].max())
    
def plot_foxsi_effarea():

    data = np.genfromtxt(data_dir + 'foxsi_effective_area.txt', skip_header = 1, delimiter = ',')
    
    energy_kev = data[:,0]
    foxsi1_cm2 = data[:,1]
    foxsi2_cm2 = data[:,2]
        
    ax = plt.subplot(111)
    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlabel('Energy [keV]')
    #ax.set_title(material.replace('_', ' ').capitalize())
    ax.set_ylabel(r'Effective Area [cm$^2$]')
		
    ax.plot(energy_kev, foxsi1_cm2, color = 'blue')
    ax.plot(energy_kev, foxsi2_cm2, color = 'red')

    ax.legend((r'FOXSI-1', r'FOXSI-2'))

    plt.show()
	
