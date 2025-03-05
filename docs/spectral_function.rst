.. _Spectral Function: 

Spectral Function
*****************

.. note:: 
    The symbols in these equations are defined on the :ref:`Definition of Symbols` page. 



The spectral radiance of :math:`\text{H}_3^+` at a particurlar wavelength can be calculated as a sum of guaussians, each representing an indivudial emission line.  

.. math::

    I(\lambda, T) = b(\lambda) + N \sum_{i=0}^{n_{lines}}\frac{I_{i}(T)}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}

where :math:`I_i(T)` is the radiance of emission line :math:`i`. It is given by: 

.. math::

    I_i(T) = \frac{ g_{ns}(2J+1)100 \times hcw_{if}A_{if}}{4\pi Q(T)}\exp{\left[-\frac{100 \times hcw_{upper}}{kT}\right]}

The :math:`\text{H}_3^+` partiation function :math:`Q(T)` is taken from Miller et al., (2013) and is expressed as a polynomial function: 

.. math::

    \log{Q} = \sum_{i=0} a_n (\log{T})^n

where the constants :math:`a_n` are provided in Table 1 of Miller et al., (2013). 

The :math:`\text{H}_3^+` line list is from Neale et al., (1996) provided in the package in a reduced form that removes emision lines with negligble radiance at temperatures < 2000 K. 