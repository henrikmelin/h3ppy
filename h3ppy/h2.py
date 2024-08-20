import os
import numpy as np

from .h3p import h3p

class h2(h3p) : 

    """
        This is an example of how to creat a `h3ppy` child class to extend the modelling and fitting functionality 
        to other molecules. Here, we use H2 line lists from Roueff et al. (2019), which is included in the package.

        In order for this to work, we need to override some of the `h3ppy` classes: 

        * The `__init__` function needs to read the correct line list. Alternatively, use the `line_list_file` keyword. 
        * `Q`: Different molecules have different partition functions, so a new one must be provided. The easiest way to do this is
          to a polynomial fit to, e.g., an ExoMol or HITRAN provided partition function. 
        * `Q_constants`: Store the polynomial coefficients here. 
        * `dQdT`: Since the partition function changes, its temperature derivative also changes. The easiest approach is to use the 
          `np.poly1D.deriv()` function. 
        * `total_emission`: The total emission calculation will also change. The simplest way to approach this is to sum the entire spectrum.

    """

    def __init__(self, line_list_file = '', **kwargs):

        # Provide the opportunity to use another line list
        if (line_list_file == '') :
            self.line_list_file = os.path.join(os.path.dirname(__file__), 'data/h2_line_list_roueff_2019_subset.txt')
        else : self.line_list_file = line_list_file

        # Set up the required parameters
        self.startup()

        # Parse any potential input
        self.parse_kwargs(kwargs)

    def Q_constants(self) :
        """
            Simple polynomial fit to the Roueff 2019 partiation function on ExoMol

            Returns:
                The :math:`H_2` partition function constants as fitted with `np.polyfit`.            
        """
        return [ 9.97306739e-10, -1.28235186e-06,  2.43332479e-02,  5.47921938e-01]
    
    def Q(self, **kwargs) :
        """
            Evaluate the partition function at a particular temperature.

            Returns:
                The evaluated :math:`H_2` partition function.
        """
        self.parse_kwargs(kwargs)

        consts = self.Q_constants()         
        Qfn    = np.poly1d(consts)

        return Qfn(self.vars['temperature'])

    def dQdT(self, **kwargs) : 	
        """
            Evaluate the temperature derivative of the partition function. 

            Returns: 
                Temperature derivative of the partition function at a particurlar temmperature. 
        """
        self.parse_kwargs(kwargs)
        consts = self.Q_constants() 

        Qfn    = np.poly1d(consts)
        dQdTfn = Qfn.deriv()

        return dQdTfn(self.vars['temperature']) 
    
    def total_emission(self, **kwargs) : 
        """
            Calculate the total amount of energy radiated by the molecule. Here, we do this by summing the whole spectrum.

            Returns: 
                The total emission of :math:`H_2`.
        """
        self.parse_kwargs(kwargs)

        # Return the wavelength integrated radiance
        return np.sum(self.calculate_line_intensities()) * self.vars['density']
        
        
        
        
        