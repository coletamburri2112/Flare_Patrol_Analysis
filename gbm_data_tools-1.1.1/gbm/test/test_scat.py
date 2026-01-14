#
#     Authors: William Cleveland (USRA),
#              Adam Goldstein (USRA) and
#              Daniel Kocevski (NASA)
#
#     Portions of the code are Copyright 2020 William Cleveland and
#     Adam Goldstein, Universities Space Research Association
#     All rights reserved.
#
#     Written for the Fermi Gamma-ray Burst Monitor (Fermi-GBM)
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import numpy as np
import os
from unittest import TestCase

from gbm.data.scat import *

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestParameter(TestCase):
    def test_attributes(self):
        param = Parameter(10.0, 0.1, 'test', units='km', support=(0, np.inf))
        self.assertEqual(param.value, 10.0)
        self.assertTupleEqual(param.uncertainty, (0.1, 0.1))
        self.assertEqual(param.name, 'test')
        self.assertEqual(param.units, 'km')
        self.assertTupleEqual(param.support, (0.0, np.inf))
        self.assertTrue(param.valid_value())
        self.assertTupleEqual(param.one_sigma_range(), (9.9, 10.1))
        self.assertTupleEqual(param.to_fits_value(), (10.0, 0.1, 0.1))
    
    def test_asymmetric(self):
        param = Parameter(10.0, (0.1, 0.5), 'test')
        self.assertTupleEqual(param.uncertainty, (0.1, 0.5))
        self.assertTupleEqual(param.one_sigma_range(), (9.9, 10.5))
        self.assertTupleEqual(param.to_fits_value(), (10.0, 0.5, 0.1))
    
    def test_invalid_value(self):
        param = Parameter(-10.0, (0.1, 0.5), 'test', support=(0.0, np.inf))
        self.assertFalse(param.valid_value())
    
    def test_errors(self):
        with self.assertRaises(ValueError):
            Parameter('hello', 0.1, 'test')
        
        with self.assertRaises(TypeError):
            Parameter(10.0, 'hello', 'test')
        
        with self.assertRaises(ValueError):
            Parameter(10.0, (1.0, 1.0, 1.0), 'test')
        
    def test_photon_flux(self):
        pflux = PhotonFlux(10.0, 0.1, (50.0, 300.0))
        self.assertTupleEqual(pflux.energy_range, (50.0, 300.0))
        self.assertEqual(pflux.name, 'Photon Flux')
        self.assertEqual(pflux.units, 'ph/cm^2/s')
        self.assertTupleEqual(pflux.support, (0.0, np.inf))

    def test_photon_fluence(self):
        pfluence = PhotonFluence(10.0, 0.1, (50.0, 300.0))
        self.assertTupleEqual(pfluence.energy_range, (50.0, 300.0))
        self.assertEqual(pfluence.name, 'Photon Fluence')
        self.assertEqual(pfluence.units, 'ph/cm^2')
        self.assertTupleEqual(pfluence.support, (0.0, np.inf))

    def test_energy_flux(self):
        eflux = EnergyFlux(1e-6, 1e-7, (50.0, 300.0))
        self.assertTupleEqual(eflux.energy_range, (50.0, 300.0))
        self.assertEqual(eflux.name, 'Energy Flux')
        self.assertEqual(eflux.units, 'erg/cm^2/s')
        self.assertTupleEqual(eflux.support, (0.0, np.inf))

    def test_energy_fluence(self):
        efluence = EnergyFluence(1e-6, 1e-7, (50.0, 300.0))
        self.assertTupleEqual(efluence.energy_range, (50.0, 300.0))
        self.assertEqual(efluence.name, 'Energy Fluence')
        self.assertEqual(efluence.units, 'erg/cm^2')
        self.assertTupleEqual(efluence.support, (0.0, np.inf))


class TestModelFit(TestCase):

    @classmethod
    def setUpClass(self):
        filename = os.path.join(data_dir, 
                                'glg_scat_all_bn170817529_flnc_comp_v02.fit')
        scat = Scat.open(filename)
        self.model_fit = scat.model_fits[0]
    
    def test_attributes(self):
        self.assertTupleEqual(self.model_fit.time_range, (-0.192, 0.064))
        self.assertEqual(self.model_fit.parameters[0].name, 'Amplitude')
        self.assertAlmostEqual(self.model_fit.parameters[0].value, 0.03199916)
        self.assertEqual(self.model_fit.parameters[1].name, 'Epeak')
        self.assertAlmostEqual(self.model_fit.parameters[1].value, 215.0943, 
                               places=4)
        self.assertEqual(self.model_fit.parameters[2].name, 'Index')
        self.assertAlmostEqual(self.model_fit.parameters[2].value, 0.1438237)
        self.assertEqual(self.model_fit.parameters[3].name, 'Pivot E =fix')
        self.assertAlmostEqual(self.model_fit.parameters[3].value, 100.0)
        self.assertAlmostEqual(self.model_fit.photon_flux.value, 2.812537, 
                               places=6)
        self.assertAlmostEqual(self.model_fit.energy_flux.value, 5.502238E-07)
        self.assertAlmostEqual(self.model_fit.photon_fluence.value, 0.7173939)
        self.assertAlmostEqual(self.model_fit.energy_fluence.value, 1.403456E-07)
        self.assertTupleEqual(self.model_fit.flux_energy_range, (10.0, 1000.0))
        self.assertEqual(self.model_fit.stat_name, 'Castor C-STAT')
        self.assertAlmostEqual(self.model_fit.stat_value, 479.6157, places=4)
        self.assertEqual(self.model_fit.dof, 479)
        self.assertAlmostEqual(self.model_fit.photon_flux_50_300.value, 
                               1.825902, places=6)
        self.assertAlmostEqual(self.model_fit.energy_fluence_50_300.value, 
                               9.852592E-08)
        self.assertAlmostEqual(self.model_fit.duration_fluence.value, 0.4657327)
        
    def test_param_list(self):
        self.assertListEqual(self.model_fit.parameter_list(), 
                             ['Amplitude', 'Epeak', 'Index', 'Pivot E =fix'])
    
    def test_to_fits(self):
        fits_row = self.model_fit.to_fits_row()
                   

class TestDetectorData(TestCase):
    
    @classmethod
    def setUpClass(self):
        filename = os.path.join(data_dir, 
                                'glg_scat_all_bn170817529_flnc_comp_v02.fit')
        scat = Scat.open(filename)
        self.det_data = scat.detectors[0]
    
    def test_attributes(self):
        self.assertEqual(self.det_data.instrument, 'GBM')
        self.assertEqual(self.det_data.detector, 'BGO_00')
        self.assertEqual(self.det_data.datatype, 'TTE')
        self.assertEqual(self.det_data.filename, 'glg_tte_b0_bn170817529_v00.fit')
        self.assertEqual(self.det_data.numchans, 128)
        self.assertTrue(self.det_data.active)
        self.assertEqual(self.det_data.response, 'glg_cspec_b0_bn170817529_v04.rsp')
        self.assertTupleEqual(self.det_data.time_range, (-0.192, 0.064))
        self.assertTupleEqual(self.det_data.energy_range, (284.65, 40108.))
        self.assertTupleEqual(self.det_data.channel_range, (3, 124))
        test_mask = np.zeros(128, dtype=bool)
        test_mask[3:125] = True
        self.assertListEqual(self.det_data.channel_mask.tolist(), test_mask.tolist())
        self.assertEqual(self.det_data.energy_edges.size, 129)
        
    def test_to_fits(self):
        fits_row = self.det_data.to_fits_row()


class TestScat(TestCase):
    
    @classmethod
    def setUpClass(self):
        filename = os.path.join(data_dir, 
                                'glg_scat_all_bn170817529_flnc_comp_v02.fit')
        self.scat = Scat.open(filename)
    
    def test_attributes(self):
        self.assertEqual(self.scat.num_detectors, 4)
        self.assertIsInstance(self.scat.detectors[0], DetectorData)
        self.assertEqual(self.scat.num_fits, 1)
        self.assertIsInstance(self.scat.model_fits[0], GbmModelFit)
        self.assertTupleEqual(tuple(self.scat.headers.keys()), 
                              ('PRIMARY', 'DETECTOR DATA', 'FIT PARAMS'))
    
    def test_add_detector(self):
        self.scat.add_detector_data(self.scat.detectors[0])
        self.scat._detectors = self.scat._detectors[:-1]

    def test_add_fit(self):
        self.scat.add_model_fit(self.scat.model_fits[0])
        self.scat._model_fits = self.scat._model_fits[:-1]

    def test_errors(self):
        with self.assertRaises(TypeError):
            self.scat.add_detector_data(1.0)
        with self.assertRaises(TypeError):
            self.scat.add_model_fit(2.0)
        with self.assertRaises(NotImplementedError):
            self.scat.write('.')