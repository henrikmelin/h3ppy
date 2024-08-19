import os
import unittest
import numpy as np

import h3ppy


class test_h3ppy(unittest.TestCase) : 

    def setUp(self) : 
        self.h3p = h3ppy.h3p()

    def test_init(self) : 
        self.assertEqual(os.path.basename(self.h3p.line_list_file), 'h3p_line_list_neale_1996_subset.txt')
        self.assertEqual(self.h3p._fit_sucess, False)
        self.assertEqual(self.h3p.vars['temperature'], 1000)

        self.assertEqual(self.h3p.nbackground, 1)
        self.assertEqual(self.h3p.noffset, 1)
        self.assertEqual(self.h3p.nbackground, 1)

    def test_model(self) : 
        self.assertAlmostEqual(np.sum(self.h3p.model()), 1.000003232973444e-17)

    def test_fit(self) : 

        # Add some noise to the default spectrum and fit it
        spec = self.h3p.add_noise(snr = 100)
        fit = self.h3p.fit(data = spec)
        vars, errs = self.h3p.get_results(verbose = False)

        self.assertEqual(len(vars), 5)
        self.assertEqual(len(errs), 5)

        self.assertAlmostEqual(vars['temperature'], 1000, places = -2)
        self.assertAlmostEqual(vars['density'], 1, places = -2)

    def test_Q(self) : 
        self.assertAlmostEqual(self.h3p.Q(T = 500), 80.58039000000002)

    def test_dQdT(self) : 
        self.assertAlmostEqual(self.h3p.dQdT(T = 500), 0.24334719375)

    def test_parse_kwargs(self) : 
        self.h3p.set(T = 400)
        self.assertEqual(self.h3p.vars['temperature'], 400)
        self.h3p.set(temperature = 600)
        self.assertEqual(self.h3p.vars['temperature'], 600)

        self.h3p.set(N = 1e12)
        self.assertEqual(self.h3p.vars['density'], 1e12)
        self.h3p.set(density = 2e13)
        self.assertEqual(self.h3p.vars['density'], 2e13)

    def test_reset_params(self) : 
        self.h3p.reset_params()
        self.assertEqual(self.h3p.vars['temperature'], 1000)
        self.assertEqual(self.h3p.vars['sigma_0'], 0)

    def test_total_emission(self) : 
        self.assertAlmostEqual(self.h3p.total_emission(T = 500, N = 1e15), 6.774289875642883e-07)

    def test_check_inputs(self) : 
        self.h3p.set(data = np.random.rand(self.h3p.wavelength.shape[0]))
        self.assertTrue(self.h3p.check_inputs())