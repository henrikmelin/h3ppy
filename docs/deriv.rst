.. _Partial Derivatives:

Partial Derivatives of the Spectral Function
********************************************

For completeness, the partial derivatives used to calculate the matricies containing the simultaneous equations described in :ref:`Approach to Fitting` are provided here.  

.. note:: 
    The symbols in these equations are defined on the :ref:`Definition of Symbols` page. 

The spectral function 
---------------------

The spectral radiance at a particurlar wavelength can be calculated as a sum of guaussians.  

.. math::

    I(\lambda, T) = b(\lambda) + N \sum_{i=0}^{n_{lines}}\frac{I_{i}(T)}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}

where

.. math::

    I(T) = \frac{ g_{ns}(2J+1)100 \times hcw_{if}A_{if}}{4\pi Q(T)}\exp{\left[-\frac{100 \times hcw_{upper}}{kT}\right]}



Partial derivative of temperature
---------------------------------
.. math:: 
    \frac{\partial I}{\partial T} = N \sum_{i=0}^{n_{lines}}\frac{1}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}\frac{\partial I_i(T)}{\partial}

where 

.. math::
    \frac{\partial I_i(T)}{\partial T} = I_i(T) \frac{100 \times hcw_{upper}}{kT^2} - I_i(T) \frac{1}{Q(T)}\frac{\partial Q(T)}{\partial T} 

where :math:`{\partial Q(T)}/{\partial T}` is given below. 


Partial derivative of partition function
----------------------------------------

.. math::
    \frac{\partial Q}{\partial T} = \mbox{log}_{10} (Q)  \times \frac{\partial}{\partial T} \left( \mbox{log}_{10} (Q) \right)

.. math::
    \frac{\partial}{\partial T} \left( \mbox{log}_{10} (Q) \right) = \sum_{n=0}^{6} na_n \left( \mbox{log}_{10} T \right)^{n-1}\frac{1}{T}

Partial derivative of column density
------------------------------------

.. math::
    \frac{\partial I}{\partial N} = \sum_{i=0}^{n_{lines}}\frac{I_{i}(T)}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}


Partial derivative of background
--------------------------------

In the case of a constant background, the partial derivative is simply: 

.. math::
    \frac{\partial I}{\partial b} = 1

However, h3ppy can also be configured to use a polynomial background function: 

.. math:: 
    b(\lambda) = \sum_i a_n \lambda^i

In this case, the partial derivative of each term becomes: 

.. math::
    \frac{\partial I}{\partial a_{i}} = \lambda^{i}

Partial dervative of the line shift
-----------------------------------

.. math:: 

    \frac{\partial I}{\partial s} = N \sum_{i=0}^{n_{lines}}\frac{I_{i}(T)}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}\frac{2(\lambda-(\lambda_i+s(\lambda)))}{2\sigma_{i}^{2}}\frac{\partial s(\lambda)}{\partial b_{i}}

Partial dervative of the line width
-----------------------------------


The line width, :math:`\sigma(\lambda)`, describes the width of the line so that the full width at half maximum (FWHM) a particular wavelength is FWHM = :math:`2\sqrt{2\log(2)}\sigma(\lambda)`. The line width derivative is:

.. math::
    \frac{\partial I}{\partial \sigma} = N \sum_{i=0}^{n_{lines}}\frac{I_{i}(T)}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}\frac{(\lambda-(\lambda_i+s(\lambda)))^2}{2\sigma_{i}^{2}} \frac{(-1)}{\sigma^3} \frac{\partial \sigma(\lambda)}{\partial a_{i}}

where

.. math:: 

    \frac{\partial \sigma(\lambda)}{\partial a_{i}} = \lambda^{i}