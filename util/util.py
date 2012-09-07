import numpy as np
import matplotlib.pyplot as plt
import sys
import platform
import datetime
from scipy import constants as con
from scipy.special import kv
import urllib2 as url
from bs4 import BeautifulSoup
from StringIO import StringIO

from scipy import interpolate

data_dir = './data/'



'''The X-ray transmission data comes from NIST 
	(http://www.nist.gov/pml/data/xraycoef/index.cfm)'''

def xray_transmission_in_air(path_length_meters, energy_kev):
	'''Provide the Xray transmission (0 to 1) in air given a path length in meters at 
	a particular energy given in keV.
	WARNING! This function is only good up to 30 keV or 30 meters. Base data does not
	go above these values.'''
	
	# the density of air at sea level
	density_cgs = 1.2754 * (1000.0/1) * (1/100.0) ** 3.0
	total_mass = density_cgs * path_length_meters * 100.0

	air_coefficients = load_mass_attenuation_coefficients(material_name='air_dry_near_sea_level')

	co = air_coefficients[1,:] * density_cgs * path_length_meters * 100.0
	
	data_energy_kev = air_coefficients[:,0]
	data_mass_coefficient = np.exp(-co)
	f = interpolate.interp1d(data_energy_kev, data_mass_coefficient)
    
	coeff = f(energy_kev)*total_mass
	transmission = np.exp(-coeff)

	return transmission

def load_mass_attenuation_coefficients(material_name='air_dry_near_sea_level'):
	'''Load the mass attenuation coefficients and mass energy-absorption coefficients
	from the data files. The allowed material names are 
	be
	air_dry_near_sea_level
	cadmium_telluride
	cesium_iodide
	gallium_arsenide
	lead_glass
	mercuric_iodide
	water_liquid
	silicon
	'''

	filename = '/Users/schriste/Dropbox/python/heroes/util/data/XrayMassCoef_' + material_name + '.txt'
	data = np.genfromtxt(filename, comments = ';', missing_values = ' ', skip_header = 8)
	
	return data
	
def plot_mass_attenuation_coefficient(material_name='air_dry_near_sea_level'):
	'''Plot the mass the mass attenuation coefficients and mass energy-absorption 
	coefficients for a named material. See load_mass_attenuation_coefficients definition
	for list of allowed materials.'''
	
	data = load_mass_attenuation_coefficients(material_name=material_name)
	
	energy_kev = data[:,0]
	mass_atten_coeff = data[:,1]
	mass_energy_atten_coeff = data[:,2]
	
	ax = plt.subplot(111)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlabel('Energy [keV]')
	ax.set_title(material_name.replace('_', ' ').capitalize())
	ax.set_ylabel(r'Mass Attenuation Coefficient [cm$^2$/g]')
		
	ax.plot(energy_kev, mass_atten_coeff)
	ax.plot(energy_kev, mass_energy_atten_coeff)

	ax.legend((r'$\mu/\rho$', r'$\mu_{en}/\rho$'))

	plt.show()

def xray_transmission(energy_kev, thickness_um, material_name='si'):
	'''Calculate the xray transmission through a material at a given energy'''
	data = load_attenuation_length(material_name = material_name)
	data_energy_kev = data[:,0]/1000.0
	data_atten_length_um = data[:,1]
	
	func_atten_length_um = interpolate.interp1d(data_energy_kev, data_atten_length_um)
	
	transmission = np.exp( -thickness_um / func_atten_length_um(energy_kev) )
	
	return transmission
	
def xray_absorption(energy_kev, thickness_um, material_name='si'):
	'''Calculate the xray absorption in a material'''
	return 1-xray_transmission(energy_kev, thickness_um, material_name=material_name)

def detector_efficiency(energy_kev, thickness_um, material_name='si'):
	'''Calculate the detector quantum efficiency (in percent) at a given energy'''
	
	return xray_absorption(energy_kev, thickness_um, material_name=material_name)*100.0

def load_attenuation_length(material_name='si'):
	filename = '/Users/schriste/Dropbox/python/heroes/util/data/' + material_name + '_xray_atten_length.txt'
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
    result = (1.e8/9.26) * gaunt_factor(energy_kev, kt) * 1/(energy_kev * np.sqrt(kt)) *n.exp(- (energy_kev / kt))

def gaunt_factor(energy_kev, kt):
    a = 0.5*energy_kev/kt
    return np.exp(a)*kv(a)



def effective_area(energy_kev):
    """Returns the HEROES effective area in cm^2 at a particular energy given in keV."""
    

def sensitivity(integration_time, de = 5, statistical_sig = 5):
    """Returns the HEROES sensitivity at a particular energy given in keV. 
    de is the width of the energy interval in keV"""
    
    energy_kev = np.arange(20,80,10)
    de = energy_kev
    det_background = np.array([2,2,2.5,3,3,3]) * 0.001
    effective_area = np.array([80,75,60,40,15,5])
    detector_efficiency = np.array([9.8, 9.2, 9.9, 9.7, 8.9, 7.7]) * 0.1
    atmo_transmission = np.array([0.26, 2.0, 3.2, 3.7, 4.2, 4.5]) * 0.1
    background_area = 8 * 0.04
    fraction_flux = 0.8
    
    a = statistical_sig ** 2 + np.sqrt(statistical_sig ** 4 + 4*statistical_sig ** 2 *
        det_background * de * background_area * integration_time)
        
    b = 2 * effective_area * de * integration_time * atmo_transmission * detector_efficiency * fraction_flux
    
    return  a/b

def msis_atmosphere_density(var, date = '2012/09/01 00:00', latitude=55, longitude=45, height=100, start=0, stop=1000, step=10):
    
    vars = [5,11] # 5 is height, 11 is density g/cm^3
    addr = 'http://omniweb.gsfc.nasa.gov/cgi/vitmo/vitmo_model.cgi'
    data = u'model=msis&year=2000&month=01&day=01&time_flag=0&hour=1.5&geo_flag=0.&latitude=55.&longitude=45.'
    data = data + u'&height=100.&profile=1&start=0.&stop=1000.&step=50.&f10_7=&f10_7_3=&ap=&format=0&'
    data = data + 'vars=0' + str(vars[0]) + '&vars=0' + str(vars[1])
    a = url.Request(addr, data)
    f = url.urlopen(a)

    print(f.read())
    f = url.urlopen(a)
    
    data = np.genfromtxt(f, skip_header = 18, skip_footer = 16)
    
    f = interpolate.interp1d(data[:,0], data[:,1])
        
    return f