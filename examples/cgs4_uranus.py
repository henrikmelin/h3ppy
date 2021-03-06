import matplotlib.pyplot as plt
import numpy as np
import h3ppy

# Read the NASA IRTF CGS$ data from Trafton et al., (1993)
data_file = 'cgs4_uranus_u01apr92.txt'
types = {'names' : ( 'w', 'i' ), 'formats' : ('f', 'f') }
dat = np.loadtxt(data_file, skiprows=4, dtype = types)

# Need to convert the instrument FOV to units of sterradian
# The pixel width is 3.1 arcsec with a slit widht of 3.1 arcsec
# Note - there are 4.2545e10 arceconds in a sterradian
spectrum = dat['i'] * 4.2545e10 / (3.1 * 3.1)
wave = dat['w']

# Make our h3ppy object :-) 
h3p = h3ppy.h3p()

# Set the wavelength and data, and use the spectral resulution to input the 
# expected line width
h3p.set(wavelength = wave, data = spectrum, R = 1300)

# We need to guess a temperature
h3p.set(temperature = 1000)

# Let h3ppy try and guess a wavelength offset
guess = h3p.guess_offset()

# Guess the density - this'll effectively scale the spectrum to the observed spectrum
# It really is not a measure of the actual density!  
guess = h3p.guess_density()

# Let h3ppy do the fitting - this will do a full five parameter fit
fit = h3p.fit(verbose = False)    

# Get the results
vars, errs = h3p.get_results()

# Plot the results! 
fig, ax = plt.subplots()
ax.set_title('First ever H$_3^+$ spectrum from Uranus: Trafton et al. (1993)')
ax.plot(wave, spectrum * 1e6, 'o', label = 'Original CGS4 Uranus spectrum')
ax.plot(wave, fit * 1e6, label = 'h3ppy fit to data')
ax.legend(frameon = False)

# Use the h3ppy helper functions for the labels
ax.set_xlabel(h3p.xlabel())
ax.set_ylabel(h3p.ylabel(prefix = '$\mu$'))
plt.savefig('../img/cgs4_uranus_fit.png')
plt.close()
