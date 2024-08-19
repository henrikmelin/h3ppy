Mathematical approach to fitting
--------------------------------

Let :math:`F(a,b,c)` be a function with three free parameters :math:`a`, :math:`b`, and :math:`c`. We want 
to fit these three parameters to a dataset :math:`f_{i}` where :math:`i` is the :math:`i^{th}` measurement. 
Whilst the three parameters may be non-linear, we can consider linear least-squares fitting with three 
simultaneous linear equations, one for each parameter. 

We define the partial derivatives of our function:

.. math:: 

    A \equiv \frac{\partial F(a,b,c)}{\partial a}

    B \equiv \frac{\partial F(a,b,c)}{\partial b}

    C \equiv \frac{\partial F(a,b,c)}{\partial c}


To determine :math:`a`, :math:`b`, and math:`c` for our dataset :math:`f`, we need to provide initial estimates for the parameters, 
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
