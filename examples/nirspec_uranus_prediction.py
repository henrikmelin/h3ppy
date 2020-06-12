import h3ppy
import matplotlib.pyplot as plt
import numpy as np

# Create the H3+ object
h3p = h3ppy.h3p()

# Define a wavelength range, e.g. typical of an observation of the H3+ Q branch
# Specify the start and end wavelengths, and the number of wavelength elements
wave = h3p.wavegen(3.94, 4.03, 1024)

# Create a H3+ model spectrum for a set of physical parameters 
# Spectral resolution R = 1200
model = h3p.model(density = 1e15, temperature = 550, R = 20000, wavelength = wave)

# Plot the model
fig, ax = plt.subplots()
ax.plot(wave, model * 1e6)
# Automagically set the labels 
ax.set_xlabel(h3p.xlabel())
ax.set_ylabel(h3p.ylabel(prefix = '\mu'))
plt.savefig('nirspec_uranus_prediction.png')
plt.close()