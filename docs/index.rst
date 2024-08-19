h3ppy
=====

h3ppy is an open source Python package for modeling and fitting :math:`H_3^+` spectra. 

Here's how to create a :math:`H_3^+` model: 

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

.. toctree::
   :maxdepth: 1
   :caption: Mathematical Description

    Spectral Function <spectral_function.rst>
    Definition of Symbols <symbols.rst>
    Approach to Fitting <math.rst>
    Partial Derivatives <deriv.rst>


.. toctree::
   :maxdepth: 2
   :caption: API Reference

    API <api.rst>
