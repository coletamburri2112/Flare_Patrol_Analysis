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
from gbm.data.drm import RSP

data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

class TestRSP(TestCase):
    filename = os.path.join(data_dir, 'glg_cspec_n4_bn120415958_v00.rsp2')
    numchans = 128
    numebins = 140
    numdrms = 12
    trigtime = 356223561.133346
    tstart1 = 3.562234321073E+08
    tstop1 = 3.562234812601E+08
    tcent1 = (tstart1+tstop1)/2.0
    tstart12 = 3.562240055561E+08
    tstop12 = 3.562240383246E+08
    tcent12 = (tstart12+tstop12)/2.0
    
    def test_attributes(self):
        rsp = RSP.open(self.filename)
        trigtime = rsp.trigtime
        self.assertAlmostEqual(trigtime, self.trigtime, places=6)
        self.assertEqual(rsp.is_gbm_file, True)
        self.assertEqual(rsp.id, '120415958')
        self.assertEqual(rsp.filename, os.path.basename(self.filename))
        self.assertEqual(rsp.is_trigger, True)
        self.assertEqual(rsp.detector, 'n4')
        self.assertEqual(rsp.datatype, 'CSPEC')
        self.assertEqual(rsp.numchans, self.numchans)
        self.assertEqual(rsp.numebins, self.numebins)
        self.assertEqual(rsp.numdrms, self.numdrms)
        self.assertAlmostEqual(rsp.tstart[0], self.tstart1-self.trigtime)
        self.assertAlmostEqual(rsp.tstop[0], self.tstop1-self.trigtime)
        self.assertAlmostEqual(rsp.tcent[0], self.tcent1-self.trigtime)
        self.assertAlmostEqual(rsp.tstart[-1], self.tstart12-self.trigtime)
        self.assertAlmostEqual(rsp.tstop[-1], self.tstop12-self.trigtime)
        self.assertAlmostEqual(rsp.tcent[-1], self.tcent12-self.trigtime)
        self.assertEqual(list(rsp.headers.keys()), 
                         ['PRIMARY', 'EBOUNDS', 'SPECRESP MATRIX1',
                          'SPECRESP MATRIX2', 'SPECRESP MATRIX3', 'SPECRESP MATRIX4',
                          'SPECRESP MATRIX5', 'SPECRESP MATRIX6', 'SPECRESP MATRIX7',
                          'SPECRESP MATRIX8', 'SPECRESP MATRIX9', 'SPECRESP MATRIX10',
                          'SPECRESP MATRIX11', 'SPECRESP MATRIX12'])
        self.assertCountEqual(rsp.photon_bin_centroids, 
                              np.sqrt(rsp.photon_bins[0]*rsp.photon_bins[1]))
        self.assertCountEqual(rsp.photon_bin_widths, 
                              rsp.photon_bins[1]-rsp.photon_bins[0])
        self.assertCountEqual(rsp.channel_centroids, 
                              np.sqrt(rsp.ebounds['E_MIN']*rsp.ebounds['E_MAX']))
        self.assertCountEqual(rsp.channel_widths, 
                              rsp.ebounds['E_MAX']-rsp.ebounds['E_MIN'])
    
    def test_extract_one(self):
        rsp = RSP.open(self.filename)
        one_drm = rsp.extract_drm(index=0)
        self.assertEqual(one_drm.numdrms, 1)
        self.assertCountEqual(one_drm.drm(0).flatten(), rsp.drm(0).flatten())
        
        one_drm = rsp.extract_drm(atime=0.0)
        self.assertEqual(one_drm.numdrms, 1)
        self.assertCountEqual(one_drm.drm(0).flatten(), rsp.drm(2).flatten())
    
    def test_interpolate(self):
        rsp = RSP.open(self.filename)
        interp_drm = rsp.interpolate(-8.0)
        for rsp_bin, interp_bin in zip(rsp.drm(2).flatten(), 
                                       interp_drm.drm(0).flatten()):
            self.assertAlmostEqual(rsp_bin, interp_bin, places=2)
        
        interp_drm = rsp.interpolate(-1e5)
        self.assertCountEqual(interp_drm.drm(0).flatten(), rsp.drm(0).flatten())
        interp_drm = rsp.interpolate(1e5)
        self.assertCountEqual(interp_drm.drm(0).flatten(), 
                              rsp.drm(rsp.numdrms-1).flatten())
    
    def test_drm_index(self):
        rsp = RSP.open(self.filename)
        idx = rsp.drm_index((0.0, 0.0))
        self.assertCountEqual(idx, [2])        
        idx = rsp.drm_index((0.0, 100.0))
        self.assertCountEqual(idx, [2,3,4])
    
    def test_nearest_drm(self):
        rsp = RSP.open(self.filename)
        drm = rsp.nearest_drm(0.0)
        self.assertCountEqual(drm.flatten(), rsp._drm_list[2].flatten())
    
    def test_photon_effarea(self):
        rsp = RSP.open(self.filename)
        effarea = rsp.photon_effective_area(index=0)
        self.assertEqual(effarea.size, 140)
        self.assertCountEqual(effarea.lo_edges, rsp.photon_bins[0])
        self.assertCountEqual(effarea.hi_edges, rsp.photon_bins[1])
        self.assertCountEqual(effarea.counts, rsp.drm(0).sum(axis=1))

    def test_channel_effarea(self):
        rsp = RSP.open(self.filename)
        effarea = rsp.channel_effective_area(index=0)
        self.assertEqual(effarea.size, 128)
        self.assertCountEqual(effarea.lo_edges, rsp.ebounds['E_MIN'])
        self.assertCountEqual(effarea.hi_edges, rsp.ebounds['E_MAX'])
        self.assertCountEqual(effarea.counts, rsp.drm(0).sum(axis=0))
   
    def test_weighted(self):
        rsp = RSP.open(self.filename)
        drm = rsp.weighted([0.5, 0.5], [-10.0, 30.0])
        test_drm = rsp._drm_list[2]*0.5 + rsp._drm_list[3]*0.5
        self.assertCountEqual(drm.flatten(), test_drm.flatten())         

    def test_write_weighted(self):
        rsp = RSP.open(self.filename)
        rsp.write_weighted(np.array([-10.0, 30.0]), np.array([30.0, 7.0]), 
                           np.array([0.5, 0.5]), data_dir, 
                           filename='glg_cspec_n4_bn120415958_v01.rsp2')
        os.remove(os.path.join(data_dir, 'glg_cspec_n4_bn120415958_v01.rsp2'))

    def test_errors(self):
        with  self.assertRaises(IOError):
            rsp = RSP.open('42.rsp')
        
        rsp = RSP.open(self.filename)
        with self.assertRaises(AssertionError):
            rsp.drm_index((100, -10.0))
        with self.assertRaises(ValueError):
            rsp.weighted([0.5, 0.5], [1,2,3.])
    
    def test_write(self):
        rsp = RSP.open(self.filename)
        one_drm = rsp.extract_drm(index=0)
        one_drm.write(data_dir, filename='test_drm.rsp')
        one_drm2 = RSP.open(os.path.join(data_dir, 'test_drm.rsp'))
        os.remove(os.path.join(data_dir, 'test_drm.rsp'))
        
        self.assertCountEqual(one_drm.drm(0).flatten(), one_drm2.drm(0).flatten())
        self.assertCountEqual(rsp.drm(0).flatten(), one_drm2.drm(0).flatten())
    
    def test_from_arrays(self):
        rsp = RSP.open(self.filename)
        minimal_rsp = RSP.from_arrays(rsp.ebounds['E_MIN'], rsp.ebounds['E_MAX'],
                                      *rsp.photon_bins, rsp.drm(0))
        rsp = RSP.open(self.filename)
        minimal_rsp = RSP.from_arrays(rsp.ebounds['E_MIN'], rsp.ebounds['E_MAX'],
                                      *rsp.photon_bins, rsp.drm(0))
        
        self.assertCountEqual(minimal_rsp.ebounds['E_MIN'], rsp.ebounds['E_MIN'])
        self.assertCountEqual(minimal_rsp.ebounds['E_MAX'], rsp.ebounds['E_MAX'])
        self.assertCountEqual(minimal_rsp.photon_bins[0], rsp.photon_bins[0])
        self.assertCountEqual(minimal_rsp.photon_bins[1], rsp.photon_bins[1])
        self.assertCountEqual(minimal_rsp.drm(0).flatten(), rsp.drm(0).flatten())
        self.assertEqual(minimal_rsp.trigtime, None)
        self.assertEqual(minimal_rsp.tstart, None)
        self.assertEqual(minimal_rsp.tstop, None)
        self.assertEqual(minimal_rsp.detector, 'all')
        
        full_rsp = RSP.from_arrays(rsp.ebounds['E_MIN'], rsp.ebounds['E_MAX'],
                                      *rsp.photon_bins, rsp.drm(0), 
                                      trigtime=rsp.trigtime, 
                                      tstart=rsp.tstart[0]+rsp.trigtime,
                                      tstop=rsp.tstop[0]+rsp.trigtime, 
                                      detnam=rsp.detector,
                                      filename='test.rsp')
        self.assertEqual(full_rsp.trigtime, rsp.trigtime)
        self.assertEqual(full_rsp.tstart, rsp.tstart[0])
        self.assertEqual(full_rsp.tstop, rsp.tstop[0])
        self.assertEqual(full_rsp.detector, rsp.detector)
    
    def test_effarea(self):
        rsp = RSP.open(self.filename)
        e1 = rsp.photon_bin_centroids[35] # ~51 keV
        e2 = rsp.photon_bin_centroids[62] # ~305 keV
        a1 = rsp.photon_effective_area(index=0).counts[35] # ~12 cm^2
        a2 = rsp.photon_effective_area(index=0).counts[62] # ~33 cm^2
        
        self.assertAlmostEqual(rsp.effective_area(e1, index=0), a1, places=0)
        self.assertAlmostEqual(rsp.effective_area(e2, index=0), a2, places=0)
        
        with self.assertRaises(TypeError):
            rsp.effective_area('hello', index=0)
        with self.assertRaises(ValueError):
            rsp.effective_area(-10.0, index=0)
        with self.assertRaises(ValueError):
            rsp.effective_area(np.array([10.0, -10.0]), index=0)

    def test_resample(self):
        rsp = RSP.open(self.filename)
        rsp_hires = rsp.resample(num_photon_bins=280)
        self.assertEqual(rsp_hires.numebins, 280)
        self.assertEqual(rsp_hires.numchans, rsp.numchans)
        
        rsp_lores = rsp.resample(num_photon_bins=70)
        self.assertEqual(rsp_lores.numebins, 70)
        self.assertEqual(rsp_lores.numchans, rsp.numchans)
        
        rsp_edges = rsp.resample(photon_bin_edges=np.array([10.0, 50.0, 300.0, 1000.0]))
        self.assertEqual(rsp_edges.numebins, 3)
        self.assertEqual(rsp_edges.numchans, rsp.numchans)
        
        with self.assertRaises(ValueError):
            rsp.resample()
        with self.assertRaises(TypeError):
            rsp.resample(num_photon_bins='hello')
        with self.assertRaises(ValueError):
            rsp.resample(num_photon_bins=-10)
        with self.assertRaises(TypeError):
            rsp.resample(photon_bin_edges=10.0)
        with self.assertRaises(ValueError):
            rsp.resample(photon_bin_edges=np.array([0.1, 10.0, 2e6]))

    def test_rebin(self):
        rsp = RSP.open(self.filename)
        
        rsp2 = rsp.rebin(factor=4)
        self.assertEqual(rsp2.numebins, 140)
        self.assertEqual(rsp2.numchans, 32)
        
        rsp3 = rsp.rebin(edge_indices=np.array([1, 10, 20, 30, 127]))
        self.assertEqual(rsp3.numebins, 140)
        self.assertEqual(rsp3.numchans, 4)
        
        with self.assertRaises(TypeError):
            rsp.rebin(factor='hello')
        with self.assertRaises(ValueError):
            rsp.rebin(factor=-4)
        with self.assertRaises(ValueError):
            rsp.rebin(factor=17)
        with self.assertRaises(TypeError):
            rsp.rebin(edge_indices=10)
        with self.assertRaises(ValueError):
            rsp.rebin(edge_indices=[1, 10, 300])            
                
