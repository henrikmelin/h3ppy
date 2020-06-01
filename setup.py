#!/usr/bin/env python

"""
Copyright (c) 2020 Henrik Melin 

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from setuptools import setup, find_packages

setup(
    name             = 'h3ppy',
    version          = '0.2',
    author           = 'Henrik Melin',
    author_email     = 'h.melin@gmail.com',
    description      = 'Model and fit H3+ spectra',
    url              = 'https://github.com/henrikmelin/h3ppy',
    keywords         = 'infrared spectroscopy H3+ modelling',
    packages         = find_packages(),
    install_requires = ['numpy'], 
    data_files       = [['data', ['data/h3p_line_list_neale_1996_subset.txt']]],
    long_description = """
# h3ppy
A python package for modelling and fitting H<sub>3</sub><sup>+</sup> spectra

## Install via pip
```
pip3 install h3ppy
```

Full [documentation on Github](https://github.com/henrikmelin/h3ppy).
    """,
    long_description_content_type="text/markdown",
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "License :: OSI Approved :: MIT License",
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
