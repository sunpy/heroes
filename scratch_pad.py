def acgaunt(wave, te_6): 
    """Calculate continuum gaunt factor using approximations of R. Mewe (18-JUN-85) to 
    full calculations of paper VI (Arnaut and Rothenflug for ion balances)"""
    
    # number of edges
    n_edg = 6

    # OVIII, OVII, NVII,  NVI ,  CVI,  CV
	# lamda - Wavelength of edges
	# temp - Max temperature of edges
	# CC - Noramlization coefficient
	# DD - exponent coefficient

    lamda = np.array([  19,   22,   25,   29,   34,   41])
    temp = np.array([2.45,  0.9,  1.7,  0.6, 1.05, 0.37])
    CC  = np.array([  3.2, 11.3, 0.69,  2.1,  3.6, 10.3])
    DD  = np.array([ 10.3,  2.5, 10.3,  2.5, 10.3,  2.5])

	gaunt_2p = np.zeros(te_6.shape, wave.shape)
	gaunt_ff = np.zeros(te_6.shape, wave.shape)
	gaunt_fb = np.zeros(te_6.shape, wave.shape)

    # Local variables for Free-bound calculation for low temperatures

    # Number of temperature regions
	n_t_l = 5	; Number of temperature regions
	# Number of wavelength regions
	n_w_l = 4	; Number of wavelength regions

;	integer*4	il,jl		! Column and row number of case table

	al = fltarr(n_t_l,n_w_l)	; al coefficients
	bl = fltarr(n_t_l,n_w_l)	; bl coefficients
	cl = fltarr(n_t_l,n_w_l)	; cl coefficients
	dl = fltarr(n_t_l,n_w_l)	; dl coefficients

	tem_liml = fltarr(n_t_l-1)	; Limits of temperature cases
	wav_liml = fltarr(n_w_l-1)	; Limits of wavelength cases

    # Normalization Coefficients
	al[*,0] = [ .248, 5.42e-12, 3.68e-4, 1.86e+3, 6.5e-3]
	al[*,1] = [ .248, 5.42e-12, 3.68e-4,    .176,   .176]
	al[*,2] = [ .248, 5.42e-12,    .323,    .323,   .323]
	al[*,3] = [.0535,    .0535,   .0535,   .0535,  .0535]

    # Exponent Coefficients
	bl[*,0] = [-1., -9.39, -4.9, -.686, -5.41]
	bl[*,1] = [-1., -9.39, -4.9,   -1.,   -1.]
	bl[*,2] = [-1., -9.39,  -1.,   -1.,   -1.]
	bl[*,3] = [-1.,   -1.,  -1.,   -1.,   -1.]

    # Exponent Coefficients
	cl[*,0] = [.158,	0.0,	0.0,	0.0,	0.0]
	cl[*,1] = [.158,	0.0,	0.0,	.233,	.233]
	cl[*,2] = [.158,	0.0,	.16,	.16,	.16]
	cl[*,3] = [.046,	.046,	.046,	.046,	.046]

    # Exponent Coefficients
	dl[*,0] = [-1.,	0.0,	0.0,	0.0,	0.0]
	dl[*,1] = [-1.,	0.0,	0.0,	-1.,	-1.]
	dl[*,2] = [-1.,	0.0,	-1.,	-1.,	-1.]
	dl[*,3] = [-1.,	-1.,	-1.,	-1.,	-1.]

    # Boundaries btwn temp. bins
	tem_liml = np.array([.015, .018, .035, .07])
	# Boundaries btwn wave  bins
	wav_liml = np.array([227.9, 504.3, 911.9])


    # Local variables for Free-bound calculation for high temperatures

    # Number of Temperature regions
	n_t = 10
	# Number of wavelength regions
	n_w = 16	

	a = np.zeros(n_t, n_w)	; a coefficients
	b = np.zeros(n_t, n_w)	; b coefficients
	c = np.zeros(n_t, n_w)	; c coefficients
	d = np.zeros(n_t, n_w)	; d coefficients

	tem_lim = np.zeros(n_t-1)	; Limits of temperature cases
	wav_lim = np.zeros(n_w-1)	; Limits of wavelength cases

    # Normalization Coefficients
	dumarr = np.zeros(10)
	a[*, 0] = [0.68, 3.73, 5.33, 14.0, 49.0, 49.0, 49.0, 4.2, 4.2,  4.2]
	a[*, 1] = [0.68, 3.73, 5.33, 14.0, 49.0, 49.0, 49.0, 4.2, 4.2, 18.4]
	a[*, 2] = [0.68, 3.73, 5.33, 14.0, 49.0, 49.0, 49.0,5.08,5.08, 18.4]
	a[*, 3] = [0.68, 3.73, 5.33, 14.0, 49.0, 49.0, 49.0,3.75,3.75,3.75]
	a[*, 4] = [0.68, 3.73, 5.33, 14.0, 49.0, 22.4, 22.4,2.12,2.12,2.12]
	a[*, 5] = [0.68, 3.73, 5.33, 14.0, 49.0, 22.4, 46.3, 46.3,6.12,6.12]
	a[*, 6] = [0.68, 3.73, 5.33, 14.0, 12.3, 12.3, 12.3, 12.3,6.12,6.12]
	a[*, 7] = [0.68, 3.73, 5.33, 14.0, 12.3, 12.3, 10.2, 10.2,4.98,4.98]
	a[*, 8] = [0.68, 3.73, 5.33, 10.2, 10.2, 10.2, 10.2, 10.2,4.98,4.98]
	a[*, 9] = [0.68, 3.73, 5.33, 10.2, 3.9,  3.9,  2.04, 2.04, 1.1,1.1]
	a[*,10] = [0.68, 3.73,.653,  1.04, 1.04, 1.04, 1.04, 1.04,1.04,1.04]
	a[*,11] = [0.68, 3.73,.653, .653,  .653, .653, 1.04, 1.04,1.04,1.04]
	a[*,12] = [0.68,3.73, .653, .653,  .653, .653, .653, .653,.653,.653]
	a[*,13] = dumarr + 0.6
	a[*,14] = dumarr + .37
	a(*,15) = dumarr + .053

    # Exponent Coefficients
	b[*, 0] = [-1.,-1.,-1.595,-.543,-1.572,-1.572,-1.572,-.82,-.82,-.82]
	b[*, 1] = [-1.,-1.,-1.595,-.543,-1.572,-1.572,-1.572,-.82,-.82,-1.33] 
	b[*, 2] = [-1.,-1.,-1.595,-.543,-1.572,-1.572,-1.572, -1., -1.,-1.33]
	b[*, 3] = [-1.,-1.,-1.595,-.543,-1.572,-1.572,-1.572, -1., -1., -1.]
	b[*, 4] = [-1.,-1.,-1.595,-.543,-1.572, -1.2, -1.2, -1., -1., -1.]
	b[*, 5] = [-1.,-1.,-1.595,-.543,-1.572,-1.2,-3.06,-3.06,-1.556,-1.556]
	b[*, 6] = [-1.,-1.,-1.595,-.543,-2.09,-2.09,-2.09,-2.09,-1.556,-1.556]
	b[*, 7] = [-1.,-1.,-1.595,-.543,-2.09,-2.09,-2.19,-2.19,-1.556,-1.556]
	b[*, 8] = [-1.,-1.,-1.595,-2.19,-2.19,-2.19,-2.19,-2.19,-1.556,-1.556]
	b[*, 9] = [-1.,-1.,-1.595,-2.19,-2.763,-2.763,-1.31,-1.31,-1.,-1.]
	b[*,10] = dumarr -1.   
	b[*,11] = dumarr -1.   
	b[*,12] = dumarr -1.   
	b[*,13] = dumarr -1.   
	b[*,14] = dumarr -1.   
	b[*,15] = dumarr -1.   

    # Exponent Coefficients
	c[*, 0] = [0.55,0.21,0.0,  0.0,-.826,-.826,-.826,   4.,4.,4.]
	c[*, 1] = [0.55,0.21,0.0,  0.0,-.826,-.826,-.826,   4.,4.,0.0]
	c[*, 2] = [0.55,0.21,0.0,  0.0,-.826,-.826,-.826,  3.9,3.9,0.0]
	c[*, 3] = [0.55,0.21,0.0,  0.0,-.826,-.826,-.826,  4.2,4.2,4.2]
	c[*, 4] = [0.55,0.21,0.0,  0.0,-.826,  0.0,  0.0,  5.6,5.6,5.6]
	c[*, 5] = [0.55,0.21,0.0,  0.0,-.826,  0.0,  0.0,  0.0,0.0,0.0]
	c[*, 6] = [0.55,0.21,0.0,  0.0,-.208,-.208,-.208,-.208,0.0,0.0]
	c[*, 7] = [0.55,0.21,0.0,  0.0,-.208,-.208,-.208,-.208,0.0,0.0]
	c[*, 8] = [0.55,0.21,0.0,-.208,-.208,-.208,-.208,-.208,0.0,0.0]
	c[*, 9] = [0.55,0.21,0.0,-.208,  0.0,  0.0,  0.0,  0.0,0.58,0.58]
	c[*,10] = [0.55,0.21,0.72,0.58, 0.58, 0.58, 0.58, 0.58,0.58,0.58]
	c[*,11] = [0.55,0.21,0.72,0.72, 0.72, 0.72, 0.58, 0.58,0.58,0.58]
	c[*,12] = [0.55,0.21,0.72,0.72, 0.72, 0.72, 0.72, 0.72,0.72,0.72]
	c[*,13] = dumarr + .55
	c[*,14] = dumarr + .158
	c[*,15] = dumarr + .05

    # Exponent Coefficients
	d[*, 0] = [-1.,-1.,0.0,0.0,-1.,-1.,-1.,-1.,-1.,-1.]
	d[*, 1] = [-1.,-1.,0.0,0.0,-1.,-1.,-1.,-1.,-1.,0.0]
	d[*, 2] = [-1.,-1.,0.0,0.0,-1.,-1.,-1.,-1.,-1.,0.0]
	d[*, 3] = [-1.,-1.,0.0,0.0,-1.,-1.,-1.,-1.,-1.,-1.]
	d[*, 4] = [-1.,-1.,0.0,0.0,-1.,0.0,0.0,-1.,-1.,-1.]
	d[*, 5] = [-1.,-1.,0.0,0.0,-1.,0.0,0.0,0.0,0.0,0.0]
	d[*, 6] = [-1.,-1.,0.0,0.0,-2.,-2.,-2.,-2.,0.0,0.0]
	d[*, 7] = [-1.,-1.,0.0,0.0,-2.,-2.,-2.,-2.,0.0,0.0]
	d[*, 8] = [-1.,-1.,0.0,-2.,-2.,-2.,-2.,-2.,0.0,0.0]
	d[*, 9] = [-1.,-1.,0.0,-2.,0.0,0.0,0.0,0.0,-1.,-1.]
	d[*,10] = dumarr - 1.
	d[*,11] = dumarr - 1.
	d[*,12] = dumarr - 1.
	d[*,13] = dumarr - 1.
	d[*,14] = dumarr - 1.
	d[*,15] = dumarr - 1.

    # Boundaries between temperature bins
	tem_lim = [.2,.258,.4,.585,1.,1.5,3.,4.5,8.] 
    # Boundaries between wavelength bins
	wav_lim = [1.4,4.6,6.1,9.1,14.2,16.8,18.6,22.5,25.3,31.6, 51.9,57.0,89.8,227.9,911.9]
	
	# Calculate Free-Free Gaunt factor		: gaunt_ff
    # (good for 1.e4 < te_6*1.e6 < 1.e9 K; 1 < wave < 1000 Ang)

    if te_6 <= 1:
        gaunt_ff = 0.29 * wave ** (0.48 * (wave ** (-0.08))) * te_6 ** (0.133*np.logo10(wave)-0.2)  
    if te_6 > 1:
	    gaunt_ff = 1.01 *wave ** (0.355 * (wave ** (-0.06))) * (te_6 / 100.) ** (0.3 * (wave ** (-0.066)))   

    # Calculate 2-photon Gaunt factor

	finit_2p = where((wave gt 19.) and (wave le 200.),fcount)

    if ((wave >= 19) + (wave <= 200)):

    	dum_2p=mkdarr(wave(finit_2p),te_6)
	    te_2p=dum_2p(1,*)
	    wv_2p=dum_2p(0,*)

	    for l in lambda: 
		    alpha = 106. / (l * (te_2p ** (-0.94))
		    
		    gaunt_2p[(wave >= 19) + (wave <= 200)] += (wv_2p gt lamda(i)) * CC(i)	$
			  * ((l[i]/wv_2p)^alpha) $
		 * sqrt(abs(cos(!pi*((lamda(i)/wv_2p)-0.5))))$
			 * ((temp(i) / te_2p)^0.45)	$
			 * 10. ** ((-DD(i)*(alog10(te_2p/temp(i)))^2)>(-37))

    # Calculate Free-Bound Gaunt factor

	low0=where(te_6 lt 0.01,lcount0)
	low1=where(((te_6 ge 0.01) and (te_6 lt 0.1)),lcount1)
	high=where(te_6 ge 0.1,hcount)

    if te_6 < 0.01: gaunt_fb = 0

    # low temperature case
    if ((te_6 >= 0.01) and (te_6 < 0.1)):

        il = indd( tem_liml,Te_6(low1))
        jl = indd( wav_liml,wave )
		ijl=mkdarr(jl,il)
		il=ijl(1,*)
		jl=ijl(0,*)

		tel = reform(rebin(te_6(low1),n_elements(te_6(low1)),$
			n_elements(wave)),1,	$
			n_elements(te_6(low1))*n_elements(wave))

    	gaunt_fb(low1,0:nwave-1) = al[il,jl] * (tel ** bl[il,jl]) * np.exp( cl[il,jl] * (tel^dl[il,jl]))

    # high temperature case Te_6 >= 0.1 M K
	if(hcount gt 0) then begin

        i = indd( tem_lim,Te_6(high))
        j = indd( wav_lim,wave )

		ij=mkdarr(j,i)
		i=ij(1,*)
		j=ij(0,*)

		teh=reform(rebin(te_6(high),n_elements(te_6(high)),$
			n_elements(wave)),1,	$
			n_elements(te_6(high))*n_elements(wave))

        gaunt_fb(high,0:nwave-1) = a[i,j] * (teh ** b[i,j]) * np.exp( c[i,j] * (teh ** d[i,j]))

	endif

	return gaunt_ff + gaunt_2p + gaunt_fb


def earth_atmosphere(height_m, DENSITY = density, TEMP = temp, PRESSURE = pressure, CGS=cgs, NUM=num, MSIS = msis

	'''Written by Steven Christe'''
;             height is expected to be in meters.
; Keywords
;             DENSITY - set density to return the density (default)
;             TEMP - to return the temperature (Celsius)
;             PRESSURE - to return the pressure (kPa)
;             NUM - set if you want number density (m^-3)
;             MSIS - set if you want realistic values above 70 km!!!!!

;gas constant in units J/kg/K
    boltzman_constant = con.k

    IF NOT keyword_set(num) THEN R = 0.286 ELSE R = 1.38d-20
    
	# the troposphere (heigh_m < 11000)
	if height_m < 11000:
		temperature = 15.04 - .00649 * height_m
        pressure = 101.29 * ((temperature + 273.1)/288.08) ** 5.256

	# the lower stratosphere (11000 < height_m < 25000)
    if height_m:
	    temperature = -56.46
        pressure = 22.65 * np.exp(1.73 - .000157 * height_m)
    
	# the upper stratosphere (height_m >25000)
    if height_m:
        temperature = -131.21 + .00299 * height_m
        pressure = 2.488 * ((temperature + 273.1)/ 216.6) ** (-11.388)
    
    density = pressure / (R * (temperature + 273))
    
    IF keyword_set(cgs) THEN density = density*(1000.)/(100.0)^3

    ;need more work here!
    restore, '~/data/msis_atmosphere_model.dat'

    h_km = result.height_km
    density_cgs = interpol(result.density_cgs, h_km*1000.0, h)

def atmospheric_absorption(height_km, look_angle = 90):
	'''Provides the atmospheric absorption of xrays from a particular height in km'''
	
    h_m = findgen(500)*1000.0
    density_cgs = earth_atmosphere(h_m,/msis)
;    PLOT, density_cgs, h_m/1d3, xtit='Density [g/cm!U3!N]', ytit='Altitude [km]', /xlog, charsize = 2.0
;    h_cm = h_m*100.0

	masscoeff = result.masscoeff
    h_cm = h_m *100.0	
    index = where(h_m GE height_km*1000.0 )
    total_mass = int_tabulated(h_cm[index], density_cgs[index], /double)
    total_mass = total_mass/sin(view_angle*!PI/180.0)
    plot_title = 'h = ' + num2str(height_km, length = 5) + ' km'
ENDIF 
