# h3ppy 😁

A python package for modelling and fitting H<sub>3</sub><sup>+</sup> spectra. Great! 

## Install via pip
```
pip3 install h3ppy
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
wave = h3p.wavegen(3.94, 4.03, 512)

# Create a H3+ model spectrum for a set of physical parameters 
model = h3p.model(density = 1e14, temperature = 1000, sigma = 0.001, wavelength = wave)

plt.plot(wave, model)
plt.xlabel('Wavelength (${\mu}m$)')
plt.ylabel('Intensity ($Wm^{-2}sr^{-1}{\mu}m^{-1}$)')
plt.savefig('example_model.png')
```
This creates the following H<sub>3</sub><sup>+</sup> spectrum: 

![Model H3+ spectra](img/example_model.png)

Neat, right?! We can now generate the spectrum for any temperature and density combination, with different wavelength coverage and at different spectral resolutions. 


## Fitting observed spectra

Here we'll simulate a pretend H<sub>3</sub><sup>+</sup> observation by adding some noise to the model spectrum above, and the we'll use `h3ppy` to fit physical parameters to it. 

```python
# Generate some noise to add to the model  
noise = np.random.normal(size = model.size) * np.max(model) / 100
pretend_data = model + noise

# Set the initial guess for the data. I'm making it different from the model input
# to show that the fit actually can converge on 'real' values
h3p.set(density = 5e13, temperature = 1300)

# Fit temperature and density to the pretend data
params_to_fit = ['temperature', 'density']
fit = h3p.fit(data = pretend_data, params_to_fit = params_to_fit)

# Get the fit variables and associated errors
vars, errs = h3p.get_results()

plt.plot(wave, pretend_data, label='Pretend data')
plt.plot(wave, fit, label = '$H_3^+$ spectrum fit')
plt.legend()
plt.savefig('example_fit.png')
```
Which produces an output in the console like:

```
[h3ppy]  Spectrum parameters:
         Temperature    = 1006.6 ± 13.0 [K]
         Column density = 9.84E+13 ± 3.05E+12 [m-2]
         ------------------------------
         sigma-0 = 1.00E-03 ± 0.00E+00
         offset-0 = 0.00E+00 ± 0.00E+00
         background-0 = 0.00E+00 ± 0.00E+00
```
Which is the same temperature and density as what we produced the model with, within the error bars of the fit. Note that since we're only fitting temperature and density, the errors on the other parameters are zero. The simulated data and fit look like:

![Model H3+ spectra](img/example_fit.png)


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
H3p = h3ppy.h3p(line_list_file = '/path/to/h3p_line_list_neale_1996_detailed.txt')
```
where the `h3p_line_list_neale_1996_detailed.txt` file is available on this GitHub directory, containing a greater number of translation lines. There is also `h3p_line_list_neale_1996_very_detailed.txt` (unzip before use). Obviously using a larger line-list will slow down all aspects of `h3ppy` (sad).

You can also concoct your own line list. The reader expects five columns in this order 

1. Angular momentum quantum number of upper level (J)
2. Wavenumber of the upper level
3. Wavenumber of the transition
4. Einstein A coefficient 
5. Spin weighting 