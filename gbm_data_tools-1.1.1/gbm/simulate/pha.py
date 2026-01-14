# pha.py: Class for simulating count spectra
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
from ..data import PHA, BAK, PHAII
from ..data.primitives import TimeEnergyBins
from .generators import *


class PhaSimulator:
    """Simulate PHA data given a modeled background spectrum, detector response,
    source spectrum, and exposure.

    Parameters:
        rsp (:class:`~gbm.data.RSP`): A detector response object
        function (:class:`~gbm.spectra.functions.Function`):
            A photon model function
        params (iterable): The parameters for the function
        exposure (float): The source exposure
        bkgd (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum
        bkgd_distrib (str): The distribution from which the background is
                            simulated; either 'Poisson' or 'Gaussian'
    """

    def __init__(self, rsp, function, params, exposure, bkgd, bkgd_distrib):
        self._rsp = rsp
        self._function = function
        self._params = params
        self._exposure = exposure
        self._src_gen = SourceSpectrumGenerator(rsp, function, params,
                                                exposure)
        self.set_background(bkgd, bkgd_distrib)

    def set_rsp(self, rsp):
        """Set/change the detector response.
        
        Args:
            rsp (:class:`~gbm.data.RSP`): A detector response object        
        """
        self._src_gen = SourceSpectrumGenerator(rsp, self._function,
                                                self._params,
                                                self._exposure)
        self._rsp = rsp

    def set_source(self, function, params, exposure):
        """Set/change the source spectrum.
        
        Args:
            function (:class:`~gbm.spectra.functions.Function`):  
                A photon model function        
            params (iterable): The parameters for the function
            exposure (float): The source exposure
        """
        self._src_gen = SourceSpectrumGenerator(self._rsp, function, params,
                                                exposure)
        self._bkgd_gen._bkgd.exposure = [self._exposure] * self._rsp.numchans
        self._function = function
        self._params = params
        self._exposure = exposure

    def set_background(self, bkgd, bkgd_distrib):
        """Set/change the background model.
        
        Args:
            bkgd (:class:`~gbm.background.BackgroundSpectrum`):
                A modeled background spectrum
            bkgd_distrib (str): The distribution from which the background is
                                simulated; either 'Poisson' or 'Gaussian'
        """
        bkgd_spectrum = BackgroundSpectrum(bkgd.rates, bkgd.rate_uncertainty,
                                           bkgd.lo_edges, bkgd.hi_edges,
                                           [self._exposure] * bkgd.size)
        if bkgd_distrib == 'Poisson':
            self._bkgd_gen = PoissonBackgroundGenerator(bkgd_spectrum)
        elif bkgd_distrib == 'Gaussian':
            self._bkgd_gen = GaussianBackgroundGenerator(bkgd_spectrum)
        else:
            raise ValueError(
                "bkgd_distrib can only be 'Poisson' or 'Gaussian'")

    def simulate_background(self, num_sims):
        """Generate simulations of the modeled background spectrum
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            list of :class:`~gbm.background.BackgroundSpectrum`:
                The deviates of the background spectrum
        """
        return [next(self._bkgd_gen) for i in range(num_sims)]

    def simulate_source(self, num_sims):
        """Generate simulations of the source spectrum
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            list of :class:`~gbm.data.primitives.EnergyBins`:
                The deviates of the source spectrum
        """
        return [next(self._src_gen) for i in range(num_sims)]

    def simulate_sum(self, num_sims):
        """Generate simulations of the background + source spectrum.
        
        Args:
            num_sims (int): Number of simulations
        
        Returns:
            list of :class:`~gbm.data.primitives.EnergyBins`:
                The deviates of the total background + source spectrum
        """
        summed = [None] * num_sims
        for i in range(num_sims):
            bkgd = next(self._bkgd_gen)
            src = next(self._src_gen)

            # since background model is formed from a rate, the background
            # "counts" won't be integers.  So we use the fractional part as
            # a probability to determine if we round up or truncate.
            bkgd_counts = bkgd.counts
            bkgd_counts[bkgd_counts < 0] = 0
            bkgd_counts_int = bkgd_counts.astype(int)
            bkgd_counts_frac = bkgd_counts - bkgd_counts_int
            extra_counts = (np.random.random(
                bkgd_counts_frac.size) > bkgd_counts_frac)
            bkgd_counts_int += extra_counts.astype(int)

            counts = bkgd_counts_int + src.counts
            summed[i] = EnergyBins(counts, src.lo_edges, src.hi_edges,
                                   src.exposure)
        return summed

    def to_bak(self, num_sims, tstart=None, tstop=None, **kwargs):
        """Produce BAK objects from simulations
        
        Args:
            num_sims (int): Number of simulations
            tstart (float, optional): The start time. If not set, then is zero.
            tstop (float, optional): Then end time. If not set, then is the exposure.
            **kwargs: Options passed to :class:`~gbm.data.BAK`
        
        Returns:
            list of :class:`~gbm.data.BAK`: The simulated BAK objects
        """
        if tstart is None:
            tstart = 0.0
        if tstop is None:
            tstop = tstart + self._exposure

        baks = self.simulate_background(num_sims)
        baks = [BAK.from_data(bak, tstart, tstop, **kwargs) for bak in baks]
        return baks

    def to_pha(self, num_sims, tstart=None, tstop=None, **kwargs):
        """Produce PHA objects of the background + source from simulations
        
        Args:
            num_sims (int): Number of simulations
            tstart (float, optional): The start time. If not set, then is zero.
            tstop (float, optional): Then end time. If not set, then is the exposure.
            **kwargs: Options passed to :class:`~gbm.data.PHA`
        
        Returns:
            list of :class:`~gbm.data.PHA`: The simulated PHA objects
        """
        if tstart is None:
            tstart = 0.0
        if tstop is None:
            tstop = tstart + self._exposure

        phas = self.simulate_sum(num_sims)
        phas = [PHA.from_data(pha, tstart, tstop, **kwargs) for pha in phas]
        return phas

    def to_phaii(self, num_sims, bin_width=None, **kwargs):
        """Produce a PHAII object by concatenating the simulations.
        
        Args:
            num_sims (int): Number of simulations
            bin_widths (float, optional): The width of each time bin.  Must be
                >= the exposure.  If not set, the is the exposure.
            **kwargs: Options passed to :class:`~gbm.data.PHAII`
        
        Returns:
            :class:`~gbm.data.PHAII`: The simulated PHAII object
        """
        if bin_width is None:
            bin_width = self._exposure
        if bin_width < self._exposure:
            raise ValueError('bin_width cannot be less than exposure')

        phas = self.simulate_sum(num_sims)
        counts = np.vstack([pha.counts for pha in phas])
        edges = np.arange(num_sims + 1) * bin_width
        data = TimeEnergyBins(counts, edges[:-1], edges[1:],
                              [self._exposure] * num_sims, phas[0].lo_edges,
                              phas[0].hi_edges)
        phaii = PHAII.from_data(data, **kwargs)
        return phaii
