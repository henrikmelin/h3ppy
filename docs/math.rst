.. _Approach to Fitting:

Mathematical Approach to Fitting
********************************

Let :math:`F(a,b,c)` be a function with three free parameters :math:`a`, :math:`b`, and :math:`c`. We want 
to fit these three parameters to a dataset :math:`f_{i}` where :math:`i` is the :math:`i^{th}` measurement. 
Whilst the three parameters may be non-linear, we can consider linear least-squares fitting with three 
simultaneous linear equations, one for each parameter. 

We define the partial derivatives of our function:

.. math:: 

    A \equiv \frac{\partial F(a,b,c)}{\partial a}

    B \equiv \frac{\partial F(a,b,c)}{\partial b}

    C \equiv \frac{\partial F(a,b,c)}{\partial c}


To determine :math:`a`, :math:`b`, and :math:`c` for our dataset :math:`f`, we need to provide initial estimates for the parameters, 
:math:`a_0`, :math:`b_0`, and :math:`c_0`. We can then calculate the difference between the observations and the initial guess at 
each spectral pixel :math:`i` contained in the observed spectrum.

.. math:: 
    \eta_i = f_{i} - F_{i}(a_{0}, b_{0}, c_{0})

We invoke Cramer's Rule (Bevington, 2003) and  define the following matrices each describing three simultaneous linear equations.

.. math:: 
    {\bf Z} = \left( \begin{array}{ccc}
    \sum A_{i}^{2} & \sum A_{i}B_{i} & \sum A_{i}C_{i} \\
    \sum A_{i}B_{i} & \sum B_{i}^{2} & \sum B_{i}C_{i}\\
    \sum A_{i}C_{i} & \sum B_{i}C_{i} & \sum C_{i}^{2}
    \end{array} \right)

.. math:: 
    {\bf A}= \left( \begin{array}{ccc}
    \sum A_{i}\eta_{i} & \sum A_{i}B_{i} & \sum A_{i}C_{i} \\
    \sum A_{i}\eta_{i} & \sum B_{i}^{2} & \sum B_{i}C_{i}\\
    \sum A_{i}\eta_{i} & \sum B_{i}C_{i} & \sum C_{i}^{2}
    \end{array} \right)

.. math:: 
    {\bf B}= \left( \begin{array}{ccc}
    \sum A_{i}^{2} & \sum A_{i}\eta_{i} & \sum A_{i}C_{i} \\
    \sum A_{i}B_{i} & \sum B_{i}\eta_{i} & \sum B_{i}C_{i}\\
    \sum A_{i}C_{i} & \sum B_{i}\eta_{i} & \sum C_{i}^{2}
    \end{array} \right)

.. math:: 
    {\bf C} = \left( \begin{array}{ccc}
    \sum A_{i}^{2} & \sum A_{i}B_{i} & \sum A_{i}\eta_{i} \\
    \sum A_{i}B_{i} & \sum B_{i}^{2} & \sum B_{i}\eta_{i}\\
    \sum A_{i}C_{i} & \sum B_{i}C_{i} & \sum C_{i}\eta_{i}
    \end{array} \right)


Now the shift from the initial guesses is given by the following expressions: 

.. math:: 
    \Delta a = \frac{ | A | }{ | Z |} 

.. math:: 
    \Delta b = \frac{ | B | }{ | Z |} 

.. math:: 
    \Delta c = \frac{ | C | }{ | Z |} 

The matrices above are re-calculated until :math:`\Delta a`, :math:`\Delta b`, and :math:`\Delta c` are sufficiently small as to change the values :math:`a`, :math:`b`, and :math:`c` insignificantly.



In the case of retrieving physical properties from an observed :math:`\text{H}_3^+` spectrum, the spectral function is given on the :ref:`Spectral Function` page. In general there are a handful of parameters that can be derived from fitting :math:`\text{H}_3^+` spectra: 

1. Temperature, :math:`T`
2. Column integrated :math:`\text{H}_3^+` density, :math:`N`
3. Some background function, :math:`b(\lambda)`
4. The width of the emission lines, :math:`\sigma(\lambda)`
5. The offset from rest wavelengths, :math:`s(\lambda)`

The last three of these can be expressed as constants or as polynomial functions of wavelength of some order. `h3ppy` will generate the matricies for any given number of parameters and perform the least squares minimisation to find the best fit. 

Uncertainties
-------------
If an observed spectrum has :math:`N` paris of wavelegth and radiance values, then the standard deviation (:math:`\mu`) of a single dependent datum is given by: 

.. math:: 
    \mu = \pm \sqrt{\frac{S^\prime}{N - \nu}}

where :math:`\nu` is the degrees of freedom of the fit and :math:`S^{\prime} = \sum_{i=0}^N{\eta_i^2}`, so that :math:`S^\prime \rightarrow S` as we iterate twoards convergence. 

The formal uncertainty (:math:`\sigma`), or standard error, of a parameter (e.g., :math:`a`) is given by: 

.. math::
    \sigma_a = \frac{\mu}{\sqrt{\frac{|Z|}{|Z_{A,A}|}}}

where :math:`|Z_{A,A}|` is the determinant of the :math:`\nu - 1` dimensional matrix formed by removing row :math:`A` and column :math:`A` from matrix :math:`Z`. This uncertainty is purely based on the pread of the data.

