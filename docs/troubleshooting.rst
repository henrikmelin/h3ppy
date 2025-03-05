.. _Troubleshooting:

Troubleshooting
***************

Things can and will go wrong. 

This page is inteded to provide some support for common issues with `h3ppy`. 

Issues with modelling a spectrum
--------------------------------

* The :math:`\text{H}_3^+` line list of Neale et al, (1996) provides wavelength covrage between about 2 to 7 :math:`\mu m`, so if you generate a model outside
  of this range it will show nothing. 

* Make sure all the required parameters are set to produce the expected results. Generally, this looks like: 

  .. code-block:: python

    # Create the h3ppy object
    h3p = h3ppy.h3p()

    # Create a wavelength range
    wave = h3p.wavegen(3.4, 3.7, 1000)

    # Set the model parameters, this is the bare minimum.
    h3p.set(T = 1000, N = 1e15, R = 2700, wave = wave)

Issues with fitting a spectrum
------------------------------
If fitting a spectrum fails, you will see an error message like this: 

.. code-block::

    ERROR - Fit failed to converge - solution is numerially unstable 

This means that `h3ppy` cannot find a stable set of parameters to describe the spectrum you have provided (see :ref:`Approach to Fitting`). 
There are a number of reasons why this happens. The most common are: 

* The starting parameters you have provided (e.g. temperature and column density) are too far from the actual solution. 
  Consider changing them to more resonable values. 
* The wavelength scale is wrong. This would mean that the :math:`\text{H}_3^+` emssion lines are not in the place that `h3ppy` expects them to be. 
  The code expects the lines to be at rest wavelength, so remove any Doppler shift from the lines. Whilst `h3ppy` will attempt to fit a wavelength 
  shift to the observed spectrum, if the difference is too large, this will unlikely work, and the fit will fail. 
* The line-width is wrong. This would mean that `h3ppy` will attempt to fit lines that are wider or narrower than the observed ones, which can cause trouble. 

It should be considered best practice to **always** generate a model before fitting to visually compare the starting model inputs for the fit (e.g. temperature, density, line-width)
to the observed spectrum. This is simply done by creating a model, before doing the fitting. 

.. code-block:: python

    # Read some data
    spectrum   = fits.getdata('some.data.fits')
    wavelength = fits.getdata('some.wavelength.fits')

    # Create the h3ppy object and provide some initial parameters
    h3p  = h3ppy.h3p()
    h3p.set(T = 500, N = 1e15, R = 2700, wave = wave, data = spectrum)

    # Have h3ppy guess the density, based on the data
    h3p.guess_density()

    # Generate a model based on the initial parameters
    model = h3p.model()

    # Visualise the initial model
    fig, ax = plt.subplots()
    ax.plot(wave, spectrum)
    ax.plot(wave, model)

    # Then do some fitting
    fit = h3p.fit()

.. note::
    The `h3p.guess_density()` method can be useful for this, it will scale the column density so that the peak of the  model 
    spectrum matches the peak of the observed spectrum. As long as the initial temperature is resonable, this can help significantly. 

Fitting a large number of spectra
---------------------------------

In general, it is common to fit large numbers of :math:`\text{H}_3^+` spectra, rather than just a single one. The internal workings of `h3ppy` will
retain the retreived parameters from the last fit, which can cause trouble if that previous fit failed. Therefore, it is commone practice to 
use the `h3p.reset_params()` to initialise all parameters before a new fit takes place. 

.. code-block:: python
    :emphasize-lines: 12

    # Read lots of some data
    spectra    = fits.getdata('some.data.fits')
    wavelength = fits.getdata('some.wavelength.fits')

    # Create the h3ppy object
    h3p  = h3ppy.h3p()

    # Iterate through all the spectra
    for spectrum in spectra : 

        # Reset the internal parameters
        h3p.reset_params()

        # Set some initial guess
        h3p.set(T = 500, N = 1e15, R = 2700, wave = wave, data = data)

        # Have h3ppy guess the density, based on the data
        h3p.guess_density()

        # Fit the spectrum
        fit = h3p.fit()

        # Rerieve the parameters
        vars, errs = h3p.get_results()

        # Check that the fit was successful 
        if (vars) : 
            ... # Store the variable 

The last lines of this code will also check if the fit was successful before storing the fit results. 



