import matplotlib.pyplot as plt
import numpy as np
import h3ppy

# Read the NASA IRTF CGS4 data from Trafton et al., (1993)
data_file = 'nirspec_jupiter.txt'
types = {'names' : ( 'w', 'i' ), 'formats' : ('f', 'f') }
dat = np.loadtxt(data_file, skiprows=4, dtype = types)
wave = dat['w']
spec = dat['i']

# Create the h3ppy object feed data into it
# The spectral resolution of NIRSPEC is ~20k
h3p = h3ppy.h3p()

# Plot the observation
title = 'Keck NIRSPEC H$_3^+$ spectrum of the southern aurora of Jupiter'
fig, ax = plt.subplots()
ax.plot(wave, spec * 1e3)
ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'), title = title)
plt.tight_layout()
plt.savefig('../img/nirspec_jupiter_data.png')
# plt.show()
plt.close()

# This function sub-divides data centered on a list of wavelengths
def subdivide(wave, spec, middles, width = 20) : 
    ret = []
    for m in middles : 
        centre = np.abs(wave - m).argmin()
        for i in range(centre - width, centre + width) : 
            ret.append(spec[i])
    return np.array(ret)
    
# The H3+ line centeres contained withing this spectral band
centers = [3.953, 3.971, 3.986, 3.9945]
cpos = np.arange(4) * 41 + 20

# Create sub-arrays, focusing on where the H3+ lines are
subspec = subdivide(wave, spec, centers)
subwave = subdivide(wave, wave, centers)

# Set the wavelength and the data
h3p.set(wavelength = subwave, data = subspec, R = 20000)

# Create a x scale for plotting 
xx      = range(len(subspec))

# Guess the density and proceed with a five parameter fit
h3p.guess_density()
fit = h3p.fit()
vars, errs = h3p.get_results()

# Plot the fit
fig, ax = plt.subplots()
ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'), xticks = cpos, title=title)
ax.set_xticklabels(centers)
ax.legend(frameon = False)
plt.tight_layout()
plt.savefig('../img/nirspec_jupiter_fit.png')
plt.close()

# See if we can improve the fit by introducing more compexity to the model
# We now add a 3rd degree polynomail for the offset, and a second degree for the 
# background and sigma 
h3p.guess_density()
h3p.set(noffset = 3, nsigma = 2, nbackground = 2)
vtf = ['temperature', 'density', 'offset-0', 'offset-1', 'offset-2', 'sigma-0',  'sigma-1', 'background-0', 'background-1']

# Here, verbose provides additiona output to the temrminal - quite useful for trouble shooting
fit2 = h3p.fit( verbose = True )
vars, errs = h3p.get_results()
    
# Plot the new (imporoved) results    
fig, ax = plt.subplots()
ax.set_title('Keck NIRSPEC H$_3^+$ spectrum of the southern aurora of Jupiter')
ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit #1')
ax.plot(xx, fit2 * 1e3, label = 'h3ppy H$_3^+$ fit #2')
ax.set_xlabel(h3p.xlabel())
ax.set_ylabel(h3p.ylabel(prefix = 'm'))
ax.set_xticks(cpos)
ax.set_xticklabels(centers)
ax.legend(frameon = False)
plt.savefig('../img/nirspec_jupiter_better.png')
plt.show()
plt.close()

wave_new = hp3.get_rest_wavelength()



exit()    

fig, ax = plt.subplots()
label = 'offset = {o0:0.5f} + ({o1:0.5f}) * wavelength'.format(o0 = vars['offset-0'], o1 = vars['offset-1'])
ax.plot(wave, h3p.offset * 1e5, label = label)
ax.set_ylabel('Wavlength offset (10$^{-5}$ $\mu$m)')
ax.set_xlabel(h3p.xlabel())
ax.legend(frameon = False)
plt.savefig('../img/nirspec_juptier_offset.png')
plt.close()

fig, ax = plt.subplots()
diff = (wave - h3p.offset) / wave
ax.plot(wave, diff)
ax.set_ylabel('Wavlength offset (10$^{-5}$ $\mu$m)')
ax.set_ylim([1 - 1e-5, 1 + 1e-5])
ax.set_xlabel(h3p.xlabel())

plt.show()