.. _Usage:

Basic Usage
***********

To create your first :math:`H_3^+` model spectrum is easy. 

.. code-block:: python

    # Create the h3ppy object
    h3p = h3ppy.h3p()

    # Generate a wavelength scale (in microns)
    wave = h3p.wavegen(3.4, 4.1, 1000)

    # Set some physical parameters
    # T - temperature
    # N - Column density
    # R - Resolving power
    # wave - wavelength scale
    h3p.set(T = 1000, N = 1e16, R = 2700, wave = wave)

    # Generate a model spectrum
    m = h3p.model()

For additional examples, check out the GitHub Pag: `https://github.com/henrikmelin/h3ppy <https://github.com/henrikmelin/h3ppy>`_.

Setting model parameters
------------------------

The `set()`, `model()`, and `fit()` methods accepts the following inputs:

* `wavelength`, `wave`, `w` - the wavelength scale on which to produce the model.  
* `data` - the observed :math:`H_3^+` spectrum
* `temperature`, `T` - the intensity of the :math:`H_3^+` spectral lines are an exponential function of the temperature. Typical ranges for the ionosphere's of the giant planets are 400 (Saturn) to 1500 K (Jupiter). 
* `density`, `N` - the column integrated :math:`H_3^+` density, this is the number of ions along the line of sight vector.
* `sigma_n` - the _n_th polynomial constant of the spectral line width (sigma)
* `offset_n` - the _n_th polynomial constant of the wavelength offset from the rest wavelength. Doppler shift and wavelength calibration errors can offset the wavelength scale. 
* `background_n` - he _n_th polynomial constant of the displacement from the zero intensity level of the spectrum
* `nsigma` - the number of polynomial constant used for the sigma.
* `noffset` - the number of polynomial constant used for the offset.
* `nbackground` - the number of polynomial constant used for the background.
