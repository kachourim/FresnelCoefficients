#!/usr/bin/python3
# -*- coding: utf-8 -*-
from numpy import *
from matplotlib.pyplot import *
import scipy.interpolate
import warnings
warnings.filterwarnings("ignore") # removes warning in case of values beyond [-1,+1] in the definition of tet_t


# General parameters
er_1	= 1					# relative permittivity of medium 1
mr_1	= 1					# relative permeability of medium 1

# Parameters of medium 2
er_2	= "Au.txt"		# relative permittivity of medium 2
mr_2	= 1					# relative permeability of medium 2

# Incidence angle(s) [°]
tet_i	= linspace(0,90,1000) 

# Wavelength(s) [m]
lam 	= linspace(200e-9,900e-9,500) 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK	= 1


# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

Z0 	=  376.730313668
teti 	= tet_i*pi/180


def Fresnel1DAng(er1, mr1, er2, mr2):

	# ================================
	# Incidence from medium 1
	# ================================
	
	# Parameters
	Z1 	= Z0*sqrt(mr1/er1)
	Z2 	= Z0*sqrt(mr2/er2)
	n1 	= sqrt(er1*mr1)
	n2	= sqrt(er2*mr2)
	tet_t	= arcsin(n1/n2*sin(teti))

	# Parallel polarization
	R_p1 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
	R_p1[isnan(R_p1)] = 1
	T_p1	= 1 - R_p1

	# Perpendicular polarization
	R_s1 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
	R_s1[isnan(R_s1)] = 1
	T_s1	= 1 - R_s1

	tet_b = arctan(n2[0]/n1[0])*180/pi
	if isreal(tet_b) and not isnan(tet_b):
		print("Incidence from medium 1: Brewster angle for Rp at {0:.2f}°".format(tet_b))

	tet_c = arcsin(n2[0]/n1[0])*180/pi
	if isreal(tet_c) and not isnan(tet_c):
		print("Incidence from medium 1: Critical angle at {0:.2f}°".format(tet_c))

	# ================================
	# Incidence from medium 2
	# ================================
	
	# Parameters
	Z1 	= Z0*sqrt(mr2/er2)
	Z2 	= Z0*sqrt(mr1/er1)
	n1 	= sqrt(er2*mr2)
	n2	= sqrt(er1*mr1)
	tet_t	= arcsin(n1/n2*sin(teti))

	# Parallel polarization
	R_p2 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
	R_p2[isnan(R_p2)] = 1
	T_p2	= 1 - R_p2

	# Perpendicular polarization
	R_s2 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
	R_s2[isnan(R_s2)] = 1
	T_s2	= 1 - R_s2

	tet_b = arctan(n2[0]/n1[0])*180/pi
	if isreal(tet_b) and not isnan(tet_b):
		print("Incidence from medium 2: Brewster angle for Rp at {0:.2f}°".format(tet_b))
		
	tet_c = arcsin(n2[0]/n1[0])*180/pi
	if isreal(tet_c) and not isnan(tet_c):
		print("Incidence from medium 2: Critical angle at {0:.2f}°".format(tet_c))
		
		
	# ================================
	# Plot
	# ================================
	fig = figure(figsize=(6,8))

	subplot(2,1,1)
	plot(tet_i,T_p1, 'r', label="Tp")
	plot(tet_i,R_p1, '--r', label="Rp")
	plot(tet_i,T_s1, 'b', label="Ts")
	plot(tet_i,R_s1, '--b', label="Rs")
	title('Incidence from medium 1')
	xlabel("Incidence angle [°]")
	ylabel("Reflection and transmission")
	legend()

	gca().set_ylim(0,1.01)
	gca().set_xlim(tet_i[0],tet_i[-1])

	subplot(2,1,2)
	plot(tet_i,T_p2, 'r', label="Tp")
	plot(tet_i,R_p2, '--r', label="Rp")
	plot(tet_i,T_s2, 'b', label="Ts")
	plot(tet_i,R_s2, '--b', label="Rs")
	title('Incidence from medium 2')
	xlabel("Incidence angle [°]")
	ylabel("Reflection and transmission")
	legend()

	gca().set_ylim(0,1.01)
	gca().set_xlim(tet_i[0],tet_i[-1])

	tight_layout()
	fig.subplots_adjust(hspace=.5)

	show()






def Fresnel1DLam(er1, mr1, er2, mr2):

	# ================================
	# Incidence from medium 1
	# ================================
	
	# Parameters
	Z1 	= Z0*sqrt(mr1/er1)
	Z2 	= Z0*sqrt(mr2/er2)
	n1 	= sqrt(er1*mr1)
	n2	= sqrt(er2*mr2)
	tet_t	= arcsin(n1/n2*sin(teti))

	# Parallel polarization
	R_p1 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
	R_p1[isnan(R_p1)] = 1
	T_p1	= 1 - R_p1

	# Perpendicular polarization
	R_s1 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
	R_s1[isnan(R_s1)] = 1
	T_s1	= 1 - R_s1

	tet_b = arctan(n2[0]/n1[0])*180/pi
	if isreal(tet_b) and not isnan(tet_b):
		print("Incidence from medium 1: Brewster angle for Rp at {0:.2f}°".format(tet_b))

	tet_c = arcsin(n2[0]/n1[0])*180/pi
	if isreal(tet_c) and not isnan(tet_c):
		print("Incidence from medium 1: Critical angle at {0:.2f}°".format(tet_c))

	# ================================
	# Incidence from medium 2
	# ================================
	
	# Parameters
	Z1 	= Z0*sqrt(mr2/er2)
	Z2 	= Z0*sqrt(mr1/er1)
	n1 	= sqrt(er2*mr2)
	n2	= sqrt(er1*mr1)
	tet_t	= arcsin(n1/n2*sin(teti))

	# Parallel polarization
	R_p2 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
	R_p2[isnan(R_p2)] = 1
	T_p2	= 1 - R_p2

	# Perpendicular polarization
	R_s2 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
	R_s2[isnan(R_s2)] = 1
	T_s2	= 1 - R_s2

	tet_b = arctan(n2[0]/n1[0])*180/pi
	if isreal(tet_b) and not isnan(tet_b):
		print("Incidence from medium 2: Brewster angle for Rp at {0:.2f}°".format(tet_b))
		
	tet_c = arcsin(n2[0]/n1[0])*180/pi
	if isreal(tet_c) and not isnan(tet_c):
		print("Incidence from medium 2: Critical angle at {0:.2f}°".format(tet_c))
		
		
	# ================================
	# Plot
	# ================================
	fig = figure(figsize=(6,8))

	subplot(2,1,1)
	plot(lam,T_p1, 'r', label="Tp")
	plot(lam,R_p1, '--r', label="Rp")
	plot(lam,T_s1, 'b', label="Ts")
	plot(lam,R_s1, '--b', label="Rs")
	title('Incidence from medium 1')
	xlabel("Wavelength [m]")
	ylabel("Reflection and transmission")
	legend()

	gca().set_ylim(0,1.01)
	gca().set_xlim(lam[0],lam[-1])

	subplot(2,1,2)
	plot(lam,T_p2, 'r', label="Tp")
	plot(lam,R_p2, '--r', label="Rp")
	plot(lam,T_s2, 'b', label="Ts")
	plot(lam,R_s2, '--b', label="Rs")
	title('Incidence from medium 2')
	xlabel("Wavelength [m]")
	ylabel("Reflection and transmission")
	legend()

	gca().set_ylim(0,1.01)
	gca().set_xlim(lam[0],lam[-1])

	tight_layout()
	fig.subplots_adjust(hspace=.5)

	show()
	
	




	
	
def Fresnel2D(Er1, Mr1, Er2, Mr2):
	
	# ================================
	# Incidence from medium 1
	# ================================

	R_p1v = []
	R_s1v = []
	
	c = 0
	for lm in lam:
		
		er1 	= Er1[c]
		mr1 	= Mr1[c]
		er2 	= Er2[c]
		mr2 	= Mr2[c]
		c 		= c + 1
		
		Z1 	= Z0*sqrt(mr1/er1)
		Z2 	= Z0*sqrt(mr2/er2)
		n1 	= sqrt(er1*mr1)
		n2	= sqrt(er2*mr2)
		tet_t	= arcsin(n1/n2*sin(teti))

		# Parallel polarization
		R_p1 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
		R_p1[isnan(R_p1)] = 1
		R_p1v.append(R_p1)
		
		# Perpendicular polarization
		R_s1 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
		R_s1[isnan(R_s1)] = 1
		R_s1v.append(R_s1)
		
	T_p1v	= 1 - array(R_p1v)
	T_s1v	= 1 - array(R_s1v)
	
	fig = figure(figsize=(10,8))
	
	subplot(2,2,1)
	imshow(R_p1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rp")
	
	subplot(2,2,2)
	imshow(T_p1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Tp")
	
	subplot(2,2,3)
	imshow(R_s1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rs")
	
	subplot(2,2,4)
	imshow(T_s1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Ts")
	
	fig.subplots_adjust(hspace=.3)
	fig.suptitle('Incidence from medium 1') 
	

	
	# ================================
	# Incidence from medium 2
	# ================================

	R_p1v = []
	R_s1v = []
	
	c = 0
	for lm in lam:
		
		er1 	= Er1[c]
		mr1 	= Mr1[c]
		er2 	= Er2[c]
		mr2 	= Mr2[c]
		c 		= c + 1
		
		Z1 	= Z0*sqrt(mr2/er2)
		Z2 	= Z0*sqrt(mr1/er1)
		n1 	= sqrt(er2*mr2)
		n2	= sqrt(er1*mr1)
		tet_t	= arcsin(n1/n2*sin(teti))

		# Parallel polarization
		R_p1 	= abs((Z2*cos(tet_t) - Z1*cos(teti))/(Z2*cos(tet_t) + Z1*cos(teti)))**2
		R_p1[isnan(R_p1)] = 1
		R_p1v.append(R_p1)
		
		# Perpendicular polarization
		R_s1 	= abs((Z2*cos(teti) - Z1*cos(tet_t))/(Z2*cos(teti) + Z1*cos(tet_t)))**2
		R_s1[isnan(R_s1)] = 1
		R_s1v.append(R_s1)
		
	T_p1v	= 1 - array(R_p1v)
	T_s1v	= 1 - array(R_s1v)
	
	fig = figure(figsize=(10,8))
	
	subplot(2,2,1)
	imshow(R_p1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rp")
	
	subplot(2,2,2)
	imshow(T_p1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Tp")
	
	subplot(2,2,3)
	imshow(R_s1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rs")
	
	subplot(2,2,4)
	imshow(T_s1v, extent = [tet_i[0], tet_i[-1], lam[0],  lam[-1]], origin='lower', vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Ts")
	
	fig.subplots_adjust(hspace=.3)
	fig.suptitle('Incidence from medium 2') 
	
	show()
	


def InterpData(d):
	
	if type(d) == str:
	
		with open(d, 'r') as f:
			data 	= loadtxt(d, comments='#')
			
			lamd 	= data[:,0]
			d1		= data[:,1]
			d2		= data[:,2]

			if NK:		# data given in terms of n and k
				rel	= (d1**2 - d2**2) + 1j*(2*d1*d2)
			else:		# data given in terms of er_r and er_i
				rel	= d1 + 1j*d2
			
		# interpolate data
		interp = scipy.interpolate.interp1d(lamd, rel)

		return interp(lam)
			
	else:
		return  full(size(lam), d)

	
	


def main():
	
	er1 	= InterpData(er_1)
	mr1 	= InterpData(mr_1)
	er2 	= InterpData(er_2)
	mr2 	= InterpData(mr_2)

	if size(lam) == 1 and size(tet_i) > 1:
		Fresnel1DAng(er1, mr1, er2, mr2) 		# if only one wavelength but several angles are specified
	elif size(lam) > 1 and size(tet_i) == 1:
		Fresnel1DLam(er1, mr1, er2, mr2) 		# if several wavelengths but only one angle are specified
	elif size(lam) > 1 and size(tet_i) > 1:
		Fresnel2D(er1, mr1, er2, mr2)	 			# if several wavelength and angles are specified
	

main()

