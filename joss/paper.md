---
title: '`h3ppy`: An open-source Python package for modelling and fitting H$_3^+$ spectra'
tags:
  - Python
  - Astronomy
  - Giant Planets
authors:
  - name: Hernik Melin
    orcid: 0000-0001-5971-2633
    affiliation: "1" 
affiliations:
 - name: Department of Engineering and Environment, Northumbria University, UK
   index: 1
date: 21 August 2024
bibliography: paper.bib

---

# Summary

`h3ppy` is an open source Python package for modelling and fitting the near-infrared spectrum of the tri-hydrogen cation H$_3^+$. This molecular ion is a main component of the charged particle ionospheres of the giant planets (Jupiter, Saturn, Uranus, and Neptune), and observations of these systems span over 30 years [@Miller:2020], using facilities such as ground-based telescopes (e.g. Keck, Very Large Telescope, NASA Infrared Telescope Facility), orbital spacecraft (e.g. Cassini and Juno), and space-based observatories (e.g., James Webb Space Telescope). By fitting the H$_3^+$ spectra, physical properties can be determined: 1) the temperature of the upper atmosphere and 2) the column integrated density of H$_3^+$ ions. The spatial and temporal distribution of these parameters reveal the processes and dynamics that govern the upper atmospheres of the giant planets, and in particular, how this region couples to both the lower atmosphere below and to the magnetic field beyond. `h3ppy` provides the tools required to both model the H$_3^+$ spectrum and perform these spectral retrievals. 

# Statement of need

`h3ppy` seeks to simplify the process of analysing H$_3^+$ spectra by providing a standardised tool for the planetary science community. It is written in Python and installation is accessible via `pip`, which makes installing the package and maintainging it very straightforward. First-time users can generate a spectrum with only a few lines of code, whilst at the same time, the code provides more advanced control of the modelling and fitting process. 

For a given temperature and H$_3^+$ column density, the code calculates the radiance of each individual H$_3^+$ transition, then distributes these in wavelength-space by giving them Gaussian line shapes, the width of which is (principally) governed by the spectral resolving power of the spectrograph. The sum of all the individual Gaussian lines then make up the model spectrum. By making use of `numpy` array features (`numpy` is the only dependency), there are significant computational gains, making the code very fast. The fitting procedure uses a hard-coded least squares approach using partial derivatives of the spectral function, invoking Cramer's Rule [@Bevington:2003]. `h3ppy` implements the H$_3^+$ line list of @Neale:1996 and the partition function of @Miller:2010. 

By making a `h3ppy` child class, the core functionality can also be used for other molecular species, and and example on how to use `h3ppy` to model and fit quadrupole molecular hydrogen spectra (H$_2$) is provided. 

The `h3ppy` code has already been used in published scientific studies, characterising the ionosphere of Uranus [@Thomas:2023] and Jupiter [@Melin:2024]. 

# Acknowledgements

HM was supported by the STFC James Webb Fellowship (ST/W001527/2) at Northumbria University, UK.

# References