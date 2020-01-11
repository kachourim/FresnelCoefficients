# Fresnel Coefficients

Computes the power reflected and transmitted from an interface between two media using Fresnel equations. 

The medium parameters may be dispersionless or dispersive. In that latter case, the parameters may be specified via a three-column txt file. The first column of the file is the wavelength (m), the second column is the real part of the relative permittivity/permeability (or the real part of the refractive index, *n*) and the third column is the imaginary part of the relative permittivity/permeability (or the imaginary part of the refractive index, *k*). See provided txt file (Au.txt) for illustration.

The different plotting options are described in the following sections, where *Rp* and *Tp* are the reflectance and transmittance for the parallel polarization (TM), and *Rs* and *Ts* are the reflectance and transmittance for the perpendicular polarization (TE). 

Note that this python script requires the `numpy`, `matplotlib` and `scipy` libraries. 

## Scattered power versus incidence angles at a single wavelength

To obtain the Fresnel coefficients for waves scattered at an interface between vacuum and glass (assuming both media are dispersionless), you may use the following configuration:

```
# Parameters of medium 1
er_1	= 1			# relative permittivity of medium 1
mr_1	= 1			# relative permeability of medium 1

# Parameters of medium 2
er_2	= 2.25			# relative permittivity of medium 2
mr_2	= 1			# relative permeability of medium 2

# Incidence angle(s) [°]
tet_i	= linspace(0,90,1000) 

# Wavelength(s)  [m]
lam 	= 1
```

Note that the exact value of the wavelength is not important in this case since the medium parameters are dispersionless. The result of this configuration is shown in the figure below. 

<img src="/images/1D_angle.png" width="1000">

The script also provides the Brewster and the critical angles. For this configuration, the script outputs:

```
Incidence from medium 1: Brewster angle for Rp at 56.31°
Incidence from medium 2: Brewster angle for Rp at 33.69°
Incidence from medium 2: Critical angle at 41.81°
```


## Scattered power versus wavelength at a single incidence angle

The script may also be used by specifying only one incidence angle and a set of wavelengths. Here is an example where the scattered power is computed at the interface between vacuum and gold, for a wavelength range between 200 nm and 900 nm, and for an incidence angle of 45°:

```
# Parameters of medium 1
er_1	= 1			# relative permittivity of medium 1
mr_1	= 1			# relative permeability of medium 1

# Parameters of medium 2
er_2	= "Au.txt"		# relative permittivity of medium 2
mr_2	= 1			# relative permeability of medium 2

# Incidence angle(s) [°]
tet_i	= 45

# Wavelength(s) [m]
lam 	=  linspace(200e-9,900e-9,500) 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK	= 1
```

Since gold is strongly dispersive within that wavelength range, the txt file containing the corresponding data is specified for *er_2*. Note that the script automatically interpolates the wavelength points not provided in the txt file. The resulting reflectance and transmittance spectra are plotted below.

<img src="/images/1D_wavelength.png" width="450">


## Scattered power versus wavelength and incidence angle

Finally, the script can also generate 2D reflectance/transmittance maps (wavelength vs incidence angle). Here is an example for an interface between vacuum and gold:

```
# Parameters of medium 1
er_1	= 1			# relative permittivity of medium 1
mr_1	= 1			# relative permeability of medium 1

# Parameters of medium 2
er_2	= "Au.txt"		# relative permittivity of medium 2
mr_2	= 1			# relative permeability of medium 2

# Incidence angle(s) [°]
tet_i	= linspace(0,90,1000) 

# Wavelength(s) [m]
lam 	= linspace(200e-9,900e-9,500) 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK	= 1
```

The corresponding result is plotted below.

<img src="/images/2D.png" width="600">

