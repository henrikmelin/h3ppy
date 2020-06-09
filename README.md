# h3ppy 😁

A python package for modelling and fitting H<sub>3</sub><sup>+</sup> spectra. Great! 

## Install via pip
```
pip3 install h3ppy
```
Or to upgrade to the latest greatest version: 
```
```
pip3 install h3ppy --upgrade
```
## Generate a model H<sub>3</sub><sup>+</sup> spectrum 


The code below generate an example spectrum: 

```python
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
model = h3p.model(density = 1e14, temperature = 1000, R = 1000, wavelength = wave)

# Plot the model
fig, ax = plt.subplots()
ax.plot(wave, model)
# Automagically set the labels 
ax.set_xlabel(h3p.xlabel())
ax.set_ylabel(h3p.ylabel())
plt.savefig('example_model.png')
plt.close() 
```
This creates the following H<sub>3</sub><sup>+</sup> spectrum: 

![Model H3+ spectra](img/example_model.png)

Neat, right?! We can now generate the spectrum for any temperature and density combination, with different wavelength coverage and at different spectral resolutions. 


## Fitting observed spectra

Here we'll simulate a pretend H<sub>3</sub><sup>+</sup> observation by adding some noise to the model spectrum above, and the we'll use `h3ppy` to fit physical parameters to it. 

```python
# Generate some noise to add to the model  
noise = np.random.normal(size = model.size) * np.max(model) / 50
pretend_data = model + noise

# Set the initial guess for the data. I'm making it different from the model input
# to show that the fit actually can converge on 'real' values
h3p.set(density = 1e13, temperature = 1300, data = pretend_data)

# Let h3ppy make an initial guess at the density
h3p.guess_density()

# Fit temperature and density to the pretend data
fit = h3p.fit()

# Get the fit variables and associated errors
vars, errs = h3p.get_results()

# Plot the model
fig, ax = plt.subplots()
ax.plot(wave, pretend_data, label='Pretend data')
ax.plot(wave, fit, label = '$H_3^+$ spectrum fit')
ax.legend()
# Automagically set the labels 
ax.set_xlabel(h3p.xlabel())
ax.set_ylabel(h3p.ylabel())
plt.savefig('../img/example_fit.png')
plt.close() 
```
Which produces an output in the console like:

```
[h3ppy]  Spectrum parameters:
         Temperature    = 1007.2 ± 8.3 [K]
         Column density = 9.81E+13 ± 1.97E+12 [m-2]
         ------------------------------
         sigma-0 = 1.59E-03 ± 6.92E-06
         offset-0 = -4.69E-07 ± 5.93E-06
         background-0 = 1.62E-08 ± 2.14E-08
```
Which is the same temperature and density as what we produced the model with, within the error bars of the fit. The fit to the simulated H<sub>3</sub><sup>+</sup> data looks like this:

![Model H3+ spectra](img/example_fit.png)


## Real world example: NASA IRTF CGS4 Uranus spectrum

There are few spectra that are as historic as this one. It's the first spectrum of H<sub>3</sub><sup>+</sup> emission from Uranus and it was taken with the NASA Infrared Telescope Facility in Hawaii in 1992, and was published by Trafton et al. (1993, Astronomical Journal, 405, 761-766). The full code for this example is contained in [`examples/cgs4_uranus.py`](examples/cgs4_uranus.py)

```python
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

```
Which produces this fit: 

![Uranus NASA IRTF H3+ spectra](img/cgs4_uranus_fit.png)

and an output in the console of:

```
[h3ppy]  Spectrum parameters:
         Temperature    = 751.5 ± 42.7 [K]
         Column density = 1.23E+15 ± 2.73E+14 [m-2]
         ------------------------------
         sigma-0 = 1.63E-03 ± 6.49E-05
         offset-0 = -9.91E-04 ± 5.64E-05
         background-0 = 2.40E-06 ± 8.77E-07
```
which is very close, within errors, that Trafton et al. (1992) got - yas! Also, note that `h3ppy` is using line lists and partition functions that weren't available in 1992! This shows that `h3ppy` can reproduce past results, which is reassuring! 


# Input parameters

The `set()`, `model()`, and `fit()` methods accepts the following inputs:


* `wavelength` - the wavelength scale on which to produce the model
* `data` - the observed H<sub>3</sub><sup>+</sup> spectrum
* `temperature` - the intensity of the H<sub>3</sub><sup>+</sup> spectral lines are an exponential function of the temperature. Typical ranges for the ionosphere's of the giant planets are 400 (Saturn) to 1500 K (Jupiter).
* `density` - the column integrated H<sub>3</sub><sup>+</sup> density, this is the number of ions along the line of sight vector.
* `sigma-n` - the _n_th polynomial constant of the spectral line width (sigma)
* `offset-n` - the _n_th polynomial constant of the wavelength offset from the rest wavelength. Doppler shift and wavelength calibration errors can offset the wavelength scale. 
* `background-n` - he _n_th polynomial constant of the displacement from the zero intensity level of the spectrum
* `nsigma` - the number of polynomial constant used for the sigma.
* `noffset` - the number of polynomial constant used for the offset.
* `nbackground` - the number of polynomial constant used for the background.


The parameters with `-n` suffix indicates that they are the nth polynomial constant. For example, if we want use the following function to describe the sigma:
```
sigma = sigma-0 + sigma-1 * wavelength + sigma-2 * wavelength^2
```
then we need to specify the following: 
```python
h3p.set(nsigma = 3, sigma-0 = 0.1, sigma-1 = 0.01, sigma-2 = 0.001) 
```

### The line width
Here we parameterise the width of the H<sub>3</sub><sup>+</sup> lines with the `sigma-n` parameter. It is related to the full width of a line profile at half maximum (FWHM) by this expression: 
```python
FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma = 2.35482 * sigma
```

### The parameters of a spectrum
The figure below illustrates how the evaluated polynomials for the offset, sigma (via the FWHM), and background at a particular wavelength determines the way that the line-intensities are distributed over across wavelength space.   

![The parameters of a spectrum](img/spectrum_parameters.png)


### Using different H<sub>3</sub><sup>+</sup> line data

If you want to use a different set of line-data you can specify a different file when you create the `h3ppy.h3p()` object. For example: 

```
h3p = h3ppy.h3p(line_list_file = '/path/to/h3p_line_list_neale_1996_detailed.txt')
```
where the `h3p_line_list_neale_1996_detailed.txt` file is available on this GitHub directory, containing a greater number of transition lines. There is also `h3p_line_list_neale_1996_very_detailed.txt` (unzip before use). Obviously using a larger line-list will slow down all aspects of `h3ppy` (sad). 

You can also concoct your own line list. The reader expects five columns in this order 

1. Angular momentum quantum number of upper level (J<sub>upper</sub>)
2. Wavenumber of the upper level (ω<sub>upper</sub> in cm<sup>-1</sup>)
3. Wavenumber of the transition (ω in cm<sup>-1</sup>)
4. Einstein A coefficient (A)
5. Spin weighting (g<sub>ns</sub>)

Note that the reader will skip the first line of the line list file. 