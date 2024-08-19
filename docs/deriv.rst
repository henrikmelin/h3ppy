Partial Derivatives of the spectral function
============================================


Partial derivative of temperature
---------------------------------
.. math:: 
    \frac{\partial I}{\partial T} = \rho \sum_{i=0}^{n_{lines}}\frac{1}{\sigma_{i}\sqrt{2\pi}}\exp{\left(-\frac{(\lambda-(\lambda_i+s(\lambda)))^{2}}{2\sigma_{i}^{2}}\right)}\frac{\partial I_i(T)}{\partial}

where 
.. math::
    \frac{\partial I_i(T)}{\partial T} = I_i(T) \frac{100 \times hcw_{upper}}{kT^2} - I_i(T) \frac{1}{Q(T)}\frac{\partial Q(T)}{\partial T} 

where :math:{\partial Q(T)}/{\partial T}` is given below. 