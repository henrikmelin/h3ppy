import os
import numpy as np
import logging

class h3p :
    """
        Setting model parameters
        ------------------------

        The `set()`, `model()`, and `fit()` methods accepts the following inputs:

        * `wavelength`, `wave`, `w` - the wavelength scale on which to produce the model.  
        * `data` - the observed :math:`H_3^+` spectrum
        * `temperature`, `T` - the intensity of the :math:`H_3^+` spectral lines are an exponential function of the temperature. Typical ranges for the ionosphere's of the giant planets are 400 (Saturn) to 1500 K (Jupiter). 
        * `density`, `N` - the column integrated :math:`H_3^+` density, this is the number of ions along the line of sight vector.
        * `sigma_n` - the nth polynomial constant of the spectral line width (sigma)
        * `offset_n` - the nth polynomial constant of the wavelength offset from the rest wavelength. Doppler shift and wavelength calibration errors can offset the wavelength scale. 
        * `background_n` - the nth polynomial constant of the displacement from the zero intensity level of the spectrum
        * `nsigma` - the number of polynomial constant used for the sigma.
        * `noffset` - the number of polynomial constant used for the offset.
        * `nbackground` - the number of polynomial constant used for the background.
        * `R` - the spectral resolving power of the instrument. This is converted in the code to a `sigma` parameter.

    """
    def __init__(self, line_list_file = '', **kwargs):
        """
            Create a `h3p` object.
        """
        # Provide the opportunity to use another line list
        if (line_list_file == '') :
            self.line_list_file = os.path.join(os.path.dirname(__file__), 'data/h3p_line_list_neale_1996_subset.txt')
        else : self.line_list_file = line_list_file

        # Set up required parameters
        self._startup()

        # Parse any potential input
        self._parse_kwargs(kwargs)

    def _startup(self) : 
        """
            Perform a number of initialisations, define constants and default limits.

        """
        # We'll work with double precision
        self.dtype = 'double'

        # Configure logging
        logging.basicConfig(format='[h3ppy] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

        #self.errors = {}

        # Physical constants
        self.k = 1.380649e-23
        self.c = 299792458.0
        self.h = 6.62607015e-34

        # Read the spectral information
        self._read_line_list()

        # This is the number of FWHM's before and after the line centre we'll
        # evaluate the intensity for. Purely for computational speedos.
        self.sigma_limit = 5

        # When to stop iterating the fitting loop. This is defined as the ratio between
        # the first parameter delta over the current.
        self.convergence_limit = 1e-2

        # The default set of variables
        self.vars = {}
        self.vars['temperature']   = 1000
        self.vars['density']       = 1.0
        self.vars['sigma_0']       = 0.001
        self.vars['offset_0']      = 0.0
        self.vars['background_0']  = 0.0
        self.wavelength = self.wavegen(3.94, 4.01, 300)
        self.data = np.array([0], dtype = self.dtype)

        # The number of terms used for these parameters
        self.nsigma        = 1
        self.noffset       = 1
        self.nbackground   = 1

        # Internal param used for speed
        self._last_temperature = 0
        self._fit_sucess = False

    def _read_line_list(self) :
        """
            Read the line data from file into a structured array. 

            For, erhm, historical reasons, the order of the line parameters must be: 

            1. Angular momentum quantum number, :math:`J`.
            2. Wavenumber of upper state, :math:`w_{upper}`.
            3. Transition frequency, :math:`w_{if}`.
            4. Einstein A coefficient, :math:`A_{if}`.
            5. Spin weightening, :math:`g_{ns}`.

            If you want to use your own line-list, specify the `line_list_file` keyword when creating the `h3ppy` object. 

        """
        types = {'names' : ( 'Ju', 'wu', 'wl', 'EA', 'gw' ), 'formats' : ('f', 'f', 'f', 'f', 'f') }
        self.line_data = np.loadtxt(self.line_list_file, skiprows=1, dtype = types)

    def _poly_fn(self, pname, nbr) :
        """
            Evaluate an polynomial function from the parameters in `self.var`.
    
            Args: 
                pname: The basename of the variable, e.g. `offset`.
                nbr: The polynomial order. 

            Returns: 
                The evaluated polynomial. 
        """
        ret = np.zeros(len(self.wavelength))
        for k in range(nbr) :
            key = pname + '_' + str(k)
            ret += self.vars[key] * np.power(self.wavelength, k)
        return ret


    def calculate_line_intensities(self, **kwargs) :
        """
            Calculate the individual line intensities per molecule at a particlar temperature. 
            See :ref:`Spectral Function` for details. 

            Returns: 
                The intensity of each individual line. 
        """
        
        self._parse_kwargs(kwargs)

        # Don't recalculate if we've just done it
        if (self._last_temperature == self.vars['temperature']) :
            return self.line_intensity

        Q = self.Q()

        exponent   = ( -1.0 * self.line_data['wu'] *  self.h *  self.c * 100.0 ) / ( self.k * self.vars['temperature'])
        intensity  = self.line_data['gw'] * (2.0 * self.line_data['Ju'] + 1)
        intensity *=  self.h *  self.c * 100.0 * self.line_data['wl']
        intensity *= np.exp(exponent) * self.line_data['EA']  / ( Q * 4.0 * np.pi )
        
        self.line_intensity = intensity
        self._last_temperature = self.vars['temperature']

        return self.line_intensity
    
    def model(self, **kwargs) :
        """
            Generate a model :math:`H_3^+` spectrum given a set of parameters. 

            .. hint:: 
                The parameters required for this function can also be set using the `set` method. 

            Args: 
                **kwargs: Physical parameters to pass to the model calculation. 
            
            Returns: 
                A spectrum of units of spectral radiance (:math:`W m^{-2} sr^{-1} \\mu m^{-1}`). 
        """

        self._parse_kwargs(kwargs)

        line_intensity = self.calculate_line_intensities()

        self.background = self._poly_fn('background', self.nbackground)

        return self._render_spectrum(line_intensity) * self.vars['density'] + self.background

    # Generate the spectrum
    def _render_spectrum(self, line_intensity, extra_fn_multiply = '', process_fn = '', **kwargs) :
        """
            Generate a model spectrum from the calculated line intensities. 

            Args: 
                line_intenisty: The array containing the individual line intensities at a particurlar tempreature. 
                process_fn: Pass an addition function to the calculation, used for evaluating the derivatives. 

            Returns:
                The spectrum of :math:`H_3^+` at a particurlar temperature and wavelength range for one molecule. 
        """
        self._parse_kwargs(kwargs)

        spectrum   = np.zeros(len(self.wavelength), dtype = self.dtype)

        self.sigma      = self._poly_fn('sigma', self.nsigma)
        self.offset     = self._poly_fn('offset', self.noffset)
        self.background = self._poly_fn('background', self.nbackground)

        # Calculate the contribution for all lines in the line-list at each wavelength
        for i, w in enumerate(self.wavelength) : 
            exponent     = -1.0 * np.power(self.wavelength[i] - ( 1e4/self.line_data['wl'] + self.offset[i]), 2) / (2.0 * np.power(self.sigma[i], 2))
            self.line_contributions = line_intensity / ( self.sigma[i] * np.sqrt(2 * np.pi) ) * np.exp(exponent)

            # The sum of all the individual contributions is the spectrun at each spectral pixel

            # We can process the spectrum by adding additional terms,
            # used mainly for the derivaties of the spectral function.
            if (process_fn == '') : spectrum[i] = np.sum(self.line_contributions)
            else : spectrum[i] = np.sum(getattr(self, process_fn)(i))

        return spectrum

    def _Q_constants(self) :
        '''
            The partition function constants from Miller et al. (2013)

            Returns: 
                The constants for a given temperature range evaluated at a certain temperature. 
        '''
        if (self.vars['temperature'] >= 100 and self.vars['temperature'] < 1800 ) :
            constants = [-1.11391, 0.0581076, 0.000302967, -2.83724e-7, 2.31119e-10, -7.15895e-14, 1.00150e-17]
        elif (self.vars['temperature'] >= 1800 and self.vars['temperature'] < 5000 ) :
            constants = [-22125.5, 51.1539, -0.0472256, 2.26131e-5, -5.85307e-9, 7.90879e-13, -4.28349e-17]
        elif (self.vars['temperature'] >= 5000 and self.vars['temperature'] < 10000 ) :
            constants = [-654293.0, 617.630, -0.237058, 4.74466e-5, -5.20566e-9, 3.05824e-13, -7.45152e-18]
        else :
            logging.error('Partition constants out of range of temperature ' + str(self.vars['temperature']) + ' K')
            return np.array([0])

        return np.array(constants, dtype = self.dtype)

    # Calculate the H3+ partition function
    def Q(self, **kwargs) :
        """
            Calculate the :math:`H_3^+` partition function, :math:`Q(T)`. 
        
            Returns:
                The evaluated partition function (Miller et al., 2013) at a particular temperature. 
        """
        self._parse_kwargs(kwargs)
        pconst = self._Q_constants()

        Q = 0.0
        for i, const in enumerate(pconst) :
            Q += const * np.power(self.vars['temperature'], np.double(i))

        return Q;

    def _dQdT(self, **kwargs) :
        """
            The temperature derivative of the partition function. See :ref:`Partial Derivatives` for more details.

            Returns: 
                The temperature derivative of the partition function.         
        """
        self._parse_kwargs(kwargs)

        vardQdT = 0.0
        for i, const in enumerate(self._Q_constants()) :
            if (i == 0) : continue
            vardQdT += float(i) * const * np.power(self.vars['temperature'], float(i - 1))
        return vardQdT

    def _dIdT(self) :
        """
            The temperature derivative of the spectral function. See :ref:`Partial Derivatives` for more details.

            Returns: 
                The temperature derivative of the spectral function.           
        """
        line_intensity = self.calculate_line_intensities()

        const1 =  self.h *  self.c * 100 / ( np.power(self.vars['temperature'], 2) *  self.k )
        Q      = self.Q()
        dQdT   = self._dQdT()

        # Caulcate the derivative for each transtion in the line-list
        dIidT  = line_intensity * const1 * self.line_data['wu']
        dIidT -= line_intensity * dQdT / Q

        return self._render_spectrum(dIidT) * self.vars['density']

    def _dIdo(self, index) :
        """
         The wavelength polynomial constants derivative of the spectral function I. See :ref:`Partial Derivatives` for more details.
        """
        line_intensity = self.calculate_line_intensities()

        self.offset_index = index
        return self._render_spectrum(line_intensity, process_fn = '_dIdo_process') * self.vars['density']

    
    def _dIdo_process(self, i) :
        """
            Alter the spectral function for `dIdo()`
        """
        numerator_deriv = 2.0 * ( self.wavelength[i] - (1e4/self.line_data['wl'] + self.offset[i]) )
        dIdc            = 1.0 * np.power(self.wavelength[i], self.offset_index)

        return self.line_contributions * numerator_deriv / (2.0 * np.power(self.sigma[i], 2)) * dIdc

    def _dIds(self, index) :
        """
            The line width derivative of the spectral function I. See :ref:`Partial Derivatives` for more details.
        """
        line_intensity = self.calculate_line_intensities()

        self.sigma_index = index
        return  self._render_spectrum(line_intensity, process_fn = '_dIds_process') * self.vars['density']

    def _dIds_process(self, i) :
        """
            Alter the spectral function for dIds()
        """

        # The exponent in the spectral function
        exponent_numerator  = -1.0 * np.power(self.wavelength[i] - ( 1e4/self.line_data['wl'] + self.offset[i]), 2)

        # The derivative of the exponent denominator
        dIds = -2.0 * exponent_numerator / (2.0 * np.power(self.sigma[i], 3))

        # The derivative of the sigma function polynmomial functions
        dsdc = np.power(self.wavelength[i], self.sigma_index)

        return ( self.line_contributions * dIds - 1.0 * self.line_contributions / self.sigma[i]) * dsdc

    def _dIdb(self, index) :
        """
        The derivative of the background polynomial constants. See :ref:`Partial Derivatives` for more details.
        """
        return np.power(self.wavelength, index)

    def _dIdN(self) :
        """
            The partial derivative of the column density. See :ref:`Partial Derivatives` for more details.
        """
        line_intensity = self.calculate_line_intensities()
        return self._render_spectrum(line_intensity)

    def set(
            self,
            temperature: float = None,
            T: float = None,
            density: float = None,
            N: float = None,
            wavelength: np.ndarray = None,
            wave: np.ndarray = None,
            w: np.ndarray = None,
            R: float = None,
            **kwargs
        ) :
        """
        Use set() to set the model parameters that you require.

        .. code-block:: python

            h3p = h3ppy.h3p()
            h3p.set(temperature = 900, density = 1e15)

        Args:
                temperature (float): The temperature (in :math:`K`).
                T (float) : Alias for temperature.
                density (float): The column integrated density (in :math:`m^{-2}`).
                N (float) : Alias for density. 
                wavelength (np.ndarray) : The wavelength scale (in :math:`\\mu m`).
                wave (np.ndarray): Alias for wavelength. 
                w (np.ndarray): Alias for wavelength.  
                R: The spectral resolution (FWHM)

        """
        parameters = {
            'temperature': temperature,
            'T': T,
            'density': density,
            'N': N,
            'wavelength': wavelength,
            'wave': wave,
            'w': wave,
            'R': R
        }

        self._parse_kwargs(parameters | kwargs)

    def _modify_vars(self, nnew, nold, var) :
        """
            Change the number of polynomial terms. For example, if you want to use a 2nd order polynomial
            to describe the background rather than a constant, then the number of terms changes, which is
            implemented here. 
        """

        # Remove polynomial terms
        if (nnew < nold) :
            for i in range(nnew + 1, nold) :
                key = var + '_' + str(i)
                del(self.vars[key])

        # Add polynomial terms, and set them to zero
        elif (nnew > nold) :
            for i in range(nold, nnew) :
                key = var + '_' + str(i)
                self.vars[key] = 0.0

    # Parse the keyword input
    def _parse_kwargs(self, kwargs) :
        """
            Parse the required model parameters. See :ref:`Usage` for details. 
        
        """
        # Prioritise processing the n variables, that deterime the number of
        # polynomial terms to use for sigma, offset, and backround.
        for key, value in kwargs.items() :

            if (key == 'nsigma') :
                self._modify_vars(value, self.nsigma, 'sigma')
                self.nsigma = value

            elif (key == 'noffset') :
                self._modify_vars(value, self.noffset, 'offset')
                self.noffset = value

            elif (key == 'nbackground') :
                self._modify_vars(value, self.nbackground, 'background')
                self.nbackground = value

        # Now set the spectrum parameters
        ok_keys = self.vars.keys()
        # ok_keys.append('nsigma', 'noffset', 'nbackground')

        for key, value in kwargs.items() :

            if (value is None) : continue

            # Allow no dash suffix for zero order polynomials
            if (key == 'sigma') : key = 'sigma_0'
            if (key == 'offset') : key = 'offset_0'
            if (key == 'background') : key = 'background_0'

            # Alternative ways to set the line-width
            if (key == 'fwhm') :
                key = 'offset_0'
                value /= np.sqrt(2 * np.pi)
            if (key == 'R') :
                key = 'sigma_0'
                value = ( np.mean(self.wavelength) / value ) / np.sqrt(2 * np.pi)

            if key in ok_keys :
                self.vars[key] = float(value)
            # Allow the shortcut T for temperature
            elif (key == 'T') : 
                self.vars['temperature'] = float(value)
            # Allow the shortcut N for density
            elif (key == 'N') : 
                self.vars['density'] = float(value)
            # Allow shorcuts wave or w for the wavelength
            elif (key == 'wavelength' or key == 'wave' or key == 'w') : 
                self.wavelength = np.array(value, dtype = self.dtype)
            elif (key == 'data') : self.data = np.array(value, dtype = self.dtype)
            elif (key in ['nsigma', 'noffset', 'nbackground']) :
                pass
            else : logging.error('Unknown set of key/values: ' + key + ' = ' + str(value))

    def fit(self, params_to_fit = '', verbose = False, niter = 14, **kwargs) :
        """
            Fit a :math:`H_3^+` spectrum.

            .. hint:: 
                The parameters required for this function can also be set using the `set` method. 

            Args: 
                params_to_fit: Specify the parameters you want to fit, e.g., if you only want to fit temperature and 
                    density `params_to_fit = ['temperature', 'density']`. By defailt it fits temperature, density, background,
                    wavelngth offset, and line width. 
                verbose: Print more detailed information. 
                niter: How may iterations of the fitting are we going to do.      

            Returns: 
                The best fit model spectrum to the data in units of spectral radiance (:math:`W m^{-2} sr^{-1} \\mu m^{-1}`). 
        """

        self._parse_kwargs(kwargs)

        # Sanity check the inputs
        if (self._check_inputs() == False) : return False

        function_map = {'temperature' : '_dIdT()', 'density' : '_dIdN()', 'offset_n' : '_dIdo(n)', 'sigma_n' : '_dIds(n)', 'background_n' : '_dIdb(n)'}

        # The default set of params to fit
        if (params_to_fit == '') :
            self.params_to_fit = self.vars.keys()
        else : self.params_to_fit = params_to_fit

        self.convergence_arrays = []

        if (verbose) : logging.info('Number of fitting iterations is set to {niter}'.format(niter = niter))

        iternbr = 0
        for i in range(niter) :

            msg = ''

            # The difference between the observation and the data
            diff = self.data - self.model()

            elements = {}
            for k, param in enumerate(self.params_to_fit) :
                if '_' in param :
                    p1, p2 = param.split('_')
                    fn = function_map[param.replace(p2, 'n')].replace('n', p2 )
                else : fn = function_map[param]
                elements[param] = np.array ( eval('self.' + fn), dtype = self.dtype)

            # Generate the Z and ABC matricies
            Z   = np.zeros([len(elements), len(elements)], dtype = self.dtype)
            ABC = np.zeros([len(elements), len(elements), len(elements)], self.dtype)
            for x, xparam in enumerate(self.params_to_fit) :
                for y, yparam in enumerate(self.params_to_fit) :
                    Z[y][x] = np.sum(elements[xparam] * elements[yparam], dtype = self.dtype)
                    for p, pparam in enumerate(self.params_to_fit) :
                        if (y == p) : diffvar = diff
                        else : diffvar =  elements[yparam]
                        ABC[p][y][x] = np.sum(elements[xparam] * diffvar, dtype = self.dtype)

            diffs = {}
            prev_vars = {}
            for p, param in enumerate(self.params_to_fit) :
                #print('Z', Z)
                #print('ABC', ABC[p])
#                print(self.dIdT())
                detABC = np.linalg.det(ABC[p])
                detZ   = np.linalg.det(Z)
                #print('detABC', detABC)
                #print('detZ', detZ)
                if (detABC == 0 or detZ == 0) :
                    msg = '|ABC| or |Z| are zero.'
                    diff = -1
                else :

                    diff =  detABC / detZ

                    prev_vars[param] = self.vars[param]
                    self.vars[param] += diff 

                    diffs[param] = diff

                # Sanity check the retrieved variables
                if (np.isfinite(diff) == False) : '|D|/|Z| returned a NaN'
                if (self.vars['temperature'] < 100) : msg = 'Temperature is less than zero'
                if (self.vars['temperature'] > 5000) : msg = 'Temperature is larger this upper boundary of h3ppy (5000 K)'
                if (self.vars['density'] < 0) : msg = 'Density is less than zero'
                self.sigma = self._poly_fn('sigma', self.nsigma)
  #              if (np.mean(self.sigma) < 0) : msg = 'Line width is negative'
                if (msg != '') :
                    self._fit_sucess = False
                    #tip = "\n        This generally happens when the line width (sigma) and/or the wavelength scale is off."
                    tip = ''
                    logging.error('Fit failed to converge - solution is numerially unstable ') # + msg + tip)
                    if (verbose) : logging.error('In this instance: {msg}'.format(msg = msg))
                    return np.full(len(self.wavelength), 0)

            # Output the intermediate values if ya want
            if (verbose) : self.print_vars()

            self._fit_sucess = True
            self.convergence_arrays.append(diffs)

            # This parameterises the level of convergence
            fracs = [ np.abs(diffs[p] / self.vars[p]) for p in self.params_to_fit ]
            converger = np.mean( fracs ) 

            if (converger < self.convergence_limit) : break
            elif (converger > 100) : break

            if (verbose) : logging.info('The converger is {con:0.2E} after {i} iterations'.format(con = converger, i = i-1))

        if (verbose) : logging.info('Fitting concluded after {i} iterations.'.format(i = i + 1))
        if (verbose) : logging.info('The final convergence condition is {lim:0.2E}.'.format(lim = converger))

        # Calculate the errors on the retrieved parameters
        diff = self.data - self.model()
        mu = np.sqrt(np.sum(np.power(diff, 2)) / (len(self.data) - len(params_to_fit)))

        # Pre populate the error dict
        self.errors = {}
        for k in self.vars.keys() :
            self.errors[k] = 0.0

        # Calculate the errors on the fitted parameters
        for p, param in enumerate(self.params_to_fit) :
            Zpp = np.delete(np.delete(Z, p, 0), p, 1)
            self.errors[param] = mu / np.sqrt(np.linalg.det(Z)  / np.linalg.det(Zpp) )

        return self.model()


    def reset_params(self) : 
        '''
        Reset all the parameters to sensible values. 
        '''
        for key in self.vars.keys() : 
            self.vars[key] = 0.0

        self.vars['temperature']   = 1000
        self.vars['density']       = 1.0

    def wavegen(self, wstart, wend, wnbr) :
        '''
        Generate a wavelength scale using start and end wavelengths and the number of wavelength elements.

        Args: 
            wstart: The wavelength at which the range should start.
            wend: The wavelength at which the range should end.
            wnbr: The number of wavelength bins.

        Returns: 
            An array containing a wavelength range.
        '''
        res = (wend - wstart) / float(wnbr)
        self.wavelength = np.arange(wstart, wend, res)
        return self.wavelength

    def total_emission(self, **kwargs) :
        '''
        Calculate the total emitted energy by H3+ as given by Miller et al., (2013)
        '''
        self._parse_kwargs(kwargs)

        if (self.vars['temperature'] >= 30 and self.vars['temperature'] < 300) :
            coeffs = [-81.9599, 0.886768, -0.0264611, 0.000462693, -4.70108e-6, 2.84979e-8,
                      -1.03090e-10, 2.13794e-13, -2.26029e-16, 8.66357e-20]
        elif (self.vars['temperature'] >= 300 and self.vars['temperature'] < 800) :
            coeffs = [-92.2048, 0.298920, -0.000962580, 1.82712e-6, -2.04420e-9, 1.24970e-12, -3.22212e-16]
        elif (self.vars['temperature'] >= 800 and self.vars['temperature'] < 1800) :
            coeffs = [-62.7016, 0.0526104, -7.22431e-5, 5.93118e-8, -2.83755e-11, 7.35415e-15, -8.01994e-19]
        elif (self.vars['temperature'] >= 1800 and self.vars['temperature'] < 5000) :
            coeffs = [-55.7672, 0.0162530, -7.68583e-6, 1.98412e-9, -2.68044e-13, 1.47026e-17]
        else :
            logging.error( str(self.vars['temperature']) + ' K is outside the allowed range 30 < T < 5000')
            return False

        logE = 0.0
        for k, coeff in enumerate(coeffs) :
            logE += coeff * np.power(self.vars['temperature'], k)
        return np.exp(logE) * self.vars['density']

    def non_LTE_scaling(self, H2_density) :
        '''
        The non-LTE scalings provided by Miller et al. (2013)

        '''
        scalings = []
        temps = [300, 800, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
        h2dens = [1e12, 1e14, 1e16, 1e18, 1e20]
        scalings[0]  = [0.0067, 0.3931, 0.9848, 0.9998, 1.0000]
        scalings[1]  = [0.0011, 0.0313, 0.7449, 0.9965, 1.0000]
        scalings[2]  = [0.0013, 0.0108, 0.4955, 0.9866, 0.9999]
        scalings[3]  = [0.0013, 0.0064, 0.3640, 0.9701, 0.9997]
        scalings[4]  = [0.0011, 0.0049, 0.2894, 0.9499, 0.9995]
        scalings[5]  = [0.0010, 0.0042, 0.2469, 0.9299, 0.9992]
        scalings[6]  = [0.0010, 0.0042, 0.2469, 0.9299, 0.9992]
        scalings[7]  = [0.0008, 0.0036, 0.2064, 0.9014, 0.9988]
        scalings[8]  = [0.0007, 0.0034, 0.1961, 0.8923, 0.9987]
        scalings[9]  = [0.0007, 0.0033, 0.1889, 0.8854, 0.9985]
        scalings[10] = [0.0007, 0.0033, 0.1836, 0.8802, 0.9985]

        return temps, h2dens, scalings

    def get_results(self, verbose = False) :
        '''
            Return the results of the spectral fit. Note that if the fit fails, the it returns `(False, False)`.

            Returns:
                vars: A dictionary containing the best fit parameters
                errs: A dictionary containing the uncertainty on the best fit parameters

        '''
        if (self._fit_sucess == False) : return False, False

        nl = "\n"
        txt  = ' Fitted parameters:' + nl
        txt += '         Temperature    = {ds:0.1f} +/- {ed:0.1f} [K]'.format(ds = self.vars['temperature'], ed = self.errors['temperature']) + nl
        txt += '         Column density = {ds:0.2E} +/-  {ed:0.2E} [m-2]'.format(ds = self.vars['density'], ed = self.errors['density']) +  nl
        txt += '         ------------------------------' + nl
        vkeys = sorted(self.vars.keys())
        for key in vkeys :
            if (key in ['temperature', 'density']) : continue
            if ( self.errors[key] == 0 ) :
                txt += '         {ea} = {ds:0.2E}'.format(ea = key, ds = self.vars[key]) + nl
            else : txt += '         {ea} = {ds:0.2E} +/- {ed:0.2E}'.format(ea = key, ds = self.vars[key], ed=self.errors[key]) +  nl

        if (verbose == True) : logging.info(txt)

        return self.vars, self.errors

    def print_vars(self) :
        '''
        Print the current set of spectrum variables held in h3ppy.
        '''
        nl = "\n"
        txt = 'Current paramter set:' + nl
        txt += '        Temperature = {ds:0.1f}'.format(ds = self.vars['temperature']) + nl
        txt += '        Density     = {ds:0.2E}'.format(ds = self.vars['density']) + nl
        vkeys = sorted(self.vars.keys())
        for k in vkeys : #, val in enumerate(sorted(self.vars)) :
            if (k in ['temperature', 'density']) : continue
            txt += '        {key} = {ds:0.2E}'.format(key = k, ds = self.vars[k]) + nl
        logging.info(txt)

        return self.vars

    def ylabel(self, label = 'Radiance', prefix = '') :
        '''
        Returns:
            The LaTeX formatted radiance ylabel for matplotlib.
        '''
        return label + ' (' + prefix + r'Wm$^{-2}$sr$^{-1}{\mu}$m$^{-1}$)'

    def xlabel(self) :
        '''
            Returns: 
                The LaTeX formatted intesity xlabel for matplotlib.
        '''
        return r'Wavelength (${\mu}$m)'

    def guess_density(self, verbose = True, **kwargs) :
        """
            Based on the peak value of the provided data, calculate the column density to provide a better
            inital guess for use in the fitting. 

            Returns:
                The modelled spectrum with the new estimated column denisty. 
        
        """
        self._parse_kwargs(kwargs)

        if (self._check_inputs() == False) : return False

        model = self.model()
        if (np.max(model) == 0) :
            logging.error('The model generated only zeros - cannot make a guess at any parameter')
            return model

        intensity_factor = (np.nanmax(self.data) - np.nanmin(self.data)) / np.nanmax(model)
        density_guess = self.vars['density'] * intensity_factor

        if (verbose) :
            logging.info('Estimated density = {ds:0.2E} m-2'.format(ds = density_guess))

        self.set(density = density_guess)

        return model * intensity_factor


    def guess_offset(self, verbose = True, **kwargs) :
        """
            Based on the peak value of the provided data and the initial guess model, calculate the wavelength shift. 

            Returns:
                The modelled spectrum with the wavelength offset applied. 
        
        """
        self._parse_kwargs(kwargs)

        if (self._check_inputs() == False) : return False

        # Generate a H3+ model and check that it is nonzero.
        model = self.model(offset = 0)
        if (np.max(model) == 0) :
            logging.error('The model generated only zeros - cannot make a guess at any parameter')
            return model

        # Calculate the wavelength diference between the pean intensity in the model
        # and the peak intensity in the data
        max_model_pos = np.argmax(model)
        max_data_pos  = np.argmax(self.data)
        offset_guess  = self.wavelength[max_data_pos] - self.wavelength[max_model_pos]

        # Check that this wavelength offset makes sense
        if (np.abs(offset_guess) > 0.5 * np.mean(self.wavelength)) :
            txt = 'The estimated offset is larger than half of the wavlength scale. This is strange. '
            txt  += 'Is the wavelength calibration that bad? Or try using a better estimate for the temperature.'
            logging.warning(txt)
            offset_guess = 0.0

        # Log the guesses
        if (verbose) :
            logging.info('Estimated offset =  {ds:0.2E} microns'.format(ds = offset_guess))

        # Feed the guesses back into the model
        self.set(offset = offset_guess)

        return self.model()

    def _check_inputs(self) :
        """
            Perform a sanity check to see if the inputs are compatible with eachother. 

            Returns: 
                A bool that indicates if things passed the test. 

        """
        # We need to have some data to guess on
        if (self.data.shape[0] < 2) :
            logging.error('ERROR - Data array is required at this stage! E.g. h3p.fit(data = uranus)')
            return False

        # Santiy check the inputs
        if (self.data.shape[0] != self.wavelength.shape[0]) :
            logging.error('ERROR - The wavelength array has a different length to the data array! ({nd} and {nw})').format(nd = len(self.data), nw = len(self.wavelength))
            return False


        return True

    def get_rest_wavelength(self) :
        '''
        Returns:     
            The offset wavelength scale where the :math:`H_3^+` lines are at rest (i.e. lab) wavelenghts.

        '''
        self.offset = self._poly_fn('offset', self.noffset)
        return self.wavelength - self.offset

    def add_noise(self, snr = 0.0, absolute = 0.0, percent = 0.0) :
        '''
        Generate a :math:`H_3^+` model with nosie added. There are three keyword parameters that can be used to describe the level of noise. The default level is 10%.

        Args: 
            snr: Signal-to-noise ratio of the noise level
            absolute: The absolute intensity value of the noise
            percent: The percentage of the maximum model value of the noise           

        Returns:     
            The modeled spectrum with the noise added to it. 
        '''    

        # Generate a model
        model = self.model()

        # Calculate the noise scaling depending on the input
        if (snr > 0) : scale = np.max(model) / snr / 2.0
        elif (absolute > 0) : scale = absolute
        elif (percent > 0) : scale = np.max(model) * percent / 100.0
        else : scale = 0.1 * np.max(model)

        # Generate the noise        
        noise = np.random.normal(size = self.wavelength.size, scale = scale).flatten()

        # Add the noise to the model
        return model + noise 
    
    
    
    
