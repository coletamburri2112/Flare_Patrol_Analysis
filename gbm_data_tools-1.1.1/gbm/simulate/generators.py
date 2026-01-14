# generators.py: Various generators for sources and backgrounds
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

from ..background.background import BackgroundSpectrum
from ..data.primitives import EnergyBins


class SimGenerator:
    """Base class for a simulation generator
    """

    def __init__(self):
        pass

    def __iter__(self):
        return self

    def __next__(self):
        return self._simulate()

    def _simulate(self):
        pass


class PoissonBackgroundGenerator(SimGenerator):
    """Simulation generator for Poisson Background.
    Once initialized, a single deviate or many deviates can be generated::
        gen = PoissonBackgroundGenerator(bkgd)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
    
    Parameters:
        bkgd (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum
    
    Yields:
        :class:`~gbm.background.BackgroundSpectrum`:
            A Poisson random deviate of the initialized spectrum
    """

    def __init__(self, bkgd):
        super().__init__()
        self._bkgd = bkgd

    def _simulate(self):
        # the poisson count deviates in each channel
        counts = np.random.poisson(self._bkgd.counts, size=(self._bkgd.size,))
        # convert to rates...
        rates = counts / self._bkgd.exposure
        rate_uncert = np.sqrt(counts) / self._bkgd.exposure
        # ...so we can populate our background spectrum
        return BackgroundSpectrum(rates, rate_uncert, self._bkgd.lo_edges,
                                  self._bkgd.hi_edges, self._bkgd.exposure)


class GaussianBackgroundGenerator(SimGenerator):
    """Simulation generator for Gaussian Background.
    Once initialized, a single deviate or many deviates can be generated::
        gen = GaussianBackgroundGenerator(bkgd)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
    
    Parameters:
        bkgd (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum

    Yields:
        :class:`~gbm.background.BackgroundSpectrum`:
            A Gaussian random deviate of the initialized spectrum
    """

    def __init__(self, bkgd):
        super().__init__()
        self._bkgd = bkgd

    def _simulate(self):
        # the gaussian rate deviates given the "centroid" rates and 
        # rate uncertainties
        counts = np.random.normal(self._bkgd.counts, self._bkgd.count_uncertainty,
                                 size=(self._bkgd.size,))
        rates = counts/self._bkgd.exposure
        return BackgroundSpectrum(rates, self._bkgd.rate_uncertainty,
                                  self._bkgd.lo_edges, self._bkgd.hi_edges,
                                  self._bkgd.exposure)


class SourceSpectrumGenerator(SimGenerator):
    """Simulation generator for a Poisson source spectrum.
    Once initialized, a single deviate or many deviates can be generated::
        gen = SourceSpectrumGenerator(rsp, function params, exposure)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
    
    Parameters:
        rsp (:class:`~gbm.data.RSP`): A detector response object
        function (:class:`~gbm.spectra.functions.Function`):
            A photon model function
        params (iterable): The parameters for the function
        exposure (float): The source exposure

    Yields:
        :class:`~gbm.data.primitives.EnergyBins`:
            A Poisson random deviate of the initialized source spectrum
    """

    def __init__(self, rsp, function, params, exposure):
        super().__init__()
        self._rates = rsp.fold_spectrum(function.fit_eval, params) * exposure
        self._rsp = rsp
        self._exposure = [exposure] * rsp.numchans

    def _simulate(self):
        counts = np.random.poisson(self._rates, size=(self._rsp.numchans,))
        return EnergyBins(counts, self._rsp.ebounds['E_MIN'],
                          self._rsp.ebounds['E_MAX'], self._exposure)


class VariablePoissonBackground(PoissonBackgroundGenerator):
    """Simulation generator for a variable Poisson Background. This
    non-homogeneous approximation allows the amplitude of the spectrum to
    be adjusted, thereby scaling the simulated counts. Once initialized, a 
    single deviate or many deviates can be generated::
        gen = VariablePoissonBackground(bkgd)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
        
        # change the amplitude to half of the initialized amplitude
        gen.amp = 0.5
    
    Parameters:
        bkgd (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum

    Attributes:
        amp (float): The amplitude, relative to initialized spectrum.  
                     Setting ``amp=1`` gives the initialized amplitude.
    
    Yields:
        :class:`~gbm.background.BackgroundSpectrum`:
            A Poisson random deviate of the spectrum
    """

    def __init__(self, bkgd):
        super().__init__(bkgd)
        self.amp = 1.0

    def _simulate(self):
        # the poisson count deviates in each channel
        counts = np.random.poisson(self._bkgd.counts * self.amp,
                                   size=(self._bkgd.size,))
        # convert to rates...
        rates = counts / self._bkgd.exposure
        rate_uncert = np.sqrt(counts) / self._bkgd.exposure
        # ...so we can populate our background spectrum
        return BackgroundSpectrum(rates, rate_uncert, self._bkgd.lo_edges,
                                  self._bkgd.hi_edges, self._bkgd.exposure)


class VariableGaussianBackground(GaussianBackgroundGenerator):
    """Simulation generator for a variable Gaussian Background. This
    non-homogeneous approximation allows the amplitude of the spectrum to
    be adjusted, thereby scaling the simulated counts. Once initialized, a 
    single deviate or many deviates can be generated::
        gen = VariableGaussianBackground(bkgd)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
        
        # change the amplitude to twice of the initialized amplitude
        gen.amp = 2.0
    
    Parameters:
        bkgd (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum

    Attributes:
        amp (float): The amplitude, relative to initialized spectrum.  
                     Setting ``amp=1`` gives the initialized amplitude.
    
    Yields:
        :class:`~gbm.background.BackgroundSpectrum`:
            A Gaussian random deviate of the spectrum
    """

    def __init__(self, bkgd):
        super().__init__(bkgd)
        self.amp = 1.0

    def _simulate(self):
        # the gaussian rate deviates given the "centroid" rates and 
        # rate uncertainties
        rates = np.random.normal(self._bkgd.rates * self.amp,
                                 self._bkgd.rate_uncertainty * self.amp,
                                 size=(self._bkgd.size,))
        rate_uncert = self._bkgd.rate_uncertainty * self.amp

        return BackgroundSpectrum(rates, rate_uncert,
                                  self._bkgd.lo_edges, self._bkgd.hi_edges,
                                  self._bkgd.exposure)


class VariableSourceSpectrumGenerator(SourceSpectrumGenerator):
    """Simulation generator for a Poisson source spectrum, efficient for
    generating deviates when the source spectrum amplitude changes.
    Once initialized, a single deviate or many deviates can be generated::
        gen = AmpFreeSourceSpectrumGenerator(rsp, function params, exposure)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
        
        # change amplitude, and generate a new deviate
        gen.amp = 0.01
        next(gen)
    
    Parameters:
        rsp (:class:`~gbm.data.RSP`): A detector response object
        function (:class:`~gbm.spectra.functions.Function`):
            A photon model function
        params (iterable): The parameters for the function
        exposure (float): The source exposure
    
    Attributes:
        amp (float): The amplitude.  This value can be set.
    
    Yields:
        :class:`~gbm.data.primitives.EnergyBins`:
            A Poisson random deviate of the initialized source spectrum
    """

    def __init__(self, rsp, function, params, exposure):
        params_temp = [1.0]
        params_temp.extend(params[1:])
        super().__init__(rsp, function, params_temp, exposure)
        self.amp = params[0]

    def _simulate(self):
        if self.amp < 0.0:
            self.amp = 0.0
        counts = np.random.poisson(self.amp * self._rates,
                                   size=(self._rsp.numchans,))
        return EnergyBins(counts, self._rsp.ebounds['E_MIN'],
                          self._rsp.ebounds['E_MAX'], self._exposure)


class EventSpectrumGenerator(SimGenerator):
    """Simulation generator producing Poisson arrival times for a source 
    spectrum during a finite slice of time. Photon losses from deadtime and
    detection/electronic processes can be accounted for by setting the
    ``min_sep > 0``.
    Once initialized, a single deviate or many deviates can be generated::
        gen = EventSpectrumGenerator(spectrum, dt)
        
        # generate a single deviate
        next(gen)
        
        # generate 10 deviates
        [next(gen) for i in range(10)]
    
    Parameters:
        count_spectrum (np.array): An array of counts in each energy channel
        dt (float): The width of the time slice in seconds
        min_sep (float, optional): The minimum possible time separation between
                                   events.  Default is 2e-6 seconds.
    
    Attributes:
        spectrum (np.array): The counts in each channel.  This value can be set.
                             The counts array will be converted to integer type.
    Yields:
        (np.array, np.array): The arrival times and energy channels for each event
    """

    def __init__(self, count_spectrum, dt, min_sep=2e-6):
        super().__init__()
        self._min_sep = min_sep
        self._dt = dt
        self._chan_nums = None
        self._beta = None
        self.spectrum = count_spectrum

    @property
    def spectrum(self):
        return self._spectrum

    @spectrum.setter
    def spectrum(self, spectrum):
        self._spectrum = spectrum.astype(int)
        if self._spectrum.sum() == 0:
            return
        
        # where do we have counts?
        chanmask = (self._spectrum > 0)
        # the 1/rate in the time slice
        self._beta = self._dt / self._spectrum.sum()

        # get the list of channels corresponding to each count
        chan_idx = np.arange(self._spectrum.size)[chanmask]
        idx = [[idx] * counts for idx, counts in
               zip(chan_idx, self._spectrum[chanmask])]
        if len(idx) > 0:
            self._chan_nums = np.concatenate(idx)
        else:
            self._chan_nums = []

    def _simulate(self):
        # no counts
        if self.spectrum.sum() == 0:
            return None

        # Simulate arrival times for each count.  Since we are simulating 
        # counts within a finite bounded window, repeat this until all arrival
        # times are within our window
        while (True):
            times = np.random.exponential(self._beta,
                                          size=(self.spectrum.sum(),))
            times = times.cumsum()
            chans = self._chan_nums
  
            # at least one event is outside our window
            if (times[-1] > self._dt):
                continue

            # If more than one event, check if all events have >= minimum spacing.
            # If there are events with spacing less than minimum spacing, we
            # have to throw away some of those events.  The reason is that we
            # are simulating the reality of recording events with real 
            # instruments and electronics, and if the event rate is high enough
            # that events are arriving faster than the detector/electronics can
            # process, we will lose some of those events.
            while (True):
                if times.size > 1:
                    diff = (times[1:] - times[:-1])
                    if (diff.min() < self._min_sep):
                        goodmask = (diff >= self._min_sep)
                        goodmask = np.concatenate(([True], goodmask))
                        times = times[goodmask]
                        chans = chans[goodmask]
                    else:
                        break
                else:
                    break

            return (times, chans)
