h3ppy
*****

`h3ppy` is an open source Python package for modelling and fitting the near-infrared spectrum of the tri-hydrogen cation :math:`\text{H}_3^+`. This molecular ion is a main component of the charged particle ionospheres of the giant planets (Jupiter, Saturn, Uranus, and Neptune), and observations of these systems span over 30 years `(e.g., Miller et al., 2020) <https://ui.adsabs.harvard.edu/abs/2020RvMP...92c5003M>`_, using facilities such as ground-based telescopes (e.g. Keck, Very Large Telescope, NASA Infrared Telescope Facility), orbital spacecraft (e.g. Cassini and Juno), and space-based observatories (e.g., James Webb Space Telescope). By fitting the :math:`\text{H}_3^+` spectra, physical properties can be determined: 1) the temperature of the upper atmosphere and 2) the column integrated density of :math:`\text{H}_3^+` ions. The spatial and temporal distribution of these parameters reveal the processes and dynamics that govern the upper atmospheres of the giant planets, and in particular, how this region couples to both the lower atmosphere below and to the magnetic field beyond. `h3ppy` provides the tools required to both model the :math:`\text{H}_3^+` spectrum and perform these spectral retrievals. 

Statement of need
-----------------

`h3ppy` seeks to simplify the process of analysing :math:`\text{H}_3^+` spectra by providing a standardised tool for the planetary science community. It is written in Python and installation is accessible via `pip`, which makes installing the package and maintainging it very straightforward. First-time users can generate a spectrum with only a few lines of code, whilst at the same time, the code provides more advanced control of the modelling and fitting process. 

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   usage
   troubleshooting
   examples

.. toctree::
   :maxdepth: 1
   :caption: Mathematical Description

   spectral_function
   symbols
   math 
   deriv

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api


Contributing to `h3ppy`
-----------------------
There are plenty of ways to contribute to `h3ppy`, and all of them are encouraged. 

* If you find a bug, please `open a new issue on GitHub <https://github.com/henrikmelin/h3ppy/issues/new>`__, providing as much detail as possible. 
* If you have great ideas for new features, `create a new issue on GitHub <https://github.com/henrikmelin/h3ppy/issues/new>`__ to suggest it. Better still,
  fork the repository and open a pull request. 

If you run into any issues using `h3ppy`, please check the :ref:`Troubleshooting` page first. You can also ask 
questions on ask a question on `GitHub <https://github.com/henrikmelin/h3ppy/discussions/new?category=q-a>`__. Lastly, please email henrik.melin@northumbria.ac.uk for any informal discussions - alwasy happy to chat. 
 


