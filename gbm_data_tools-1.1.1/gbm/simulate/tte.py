# tte.py: Class for simulating TTE (event) data
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
import warnings
from ..data import TTE
from ..data.primitives import EventList
from .generators import *


class TteSourceSimulator:
    """Simulate TTE or EventList data for a source spectrum given a detector 
    response, spectral model and time profile model.  The spectral shape is 
    fixed throughout the time profile of the signal, but the amplitude of the 
    spectrum is time-dependent, set by the time profile function.

    Parameters:
        rsp (:class:`~gbm.data.RSP`): A detector response object
        spec_func (:class:`~gbm.spectra.functions.Function`):
            A photon model function
        spec_params (iterable): The parameters for the function
        time_func (<function>): A time profile function
        time_params (iterable): Parameters for the time profile function
        sample_period (float, optional): The sampling period of the simulator
            in seconds. Default is 0.001. The simulator will produce arrival 
            times consistent with a spectrum over a finite time slice.  This 
            time slice should be short enough to approximate a non-homogeneous 
            Poisson process, but long enough to allow for a tractable 
            computation time.
        deadtime (float, optional): The dead time in seconds for each recorded 
                                    count during which another count cannot be 
                                    recorded. Default is 2e-6 s.
    """

    def __init__(self, rsp, spec_func, spec_params, time_func, time_params,
                 sample_period=0.001, deadtime=2e-6):
        
        warnings.filterwarnings("ignore", category=UserWarning)
        
        if sample_period <= 0.0:
            raise ValueError('Sample period must be positive')
        if deadtime < 0.0:
            raise ValueError('Deadtime must be non-negative')
        
        self._rsp = rsp
        self._spec_func = spec_func
        self._spec_params = spec_params
        self._time_func = time_func
        self._time_params = time_params
        self._sample_per = sample_period
        self._spec_gen = VariableSourceSpectrumGenerator(rsp, spec_func,
                                                         spec_params,
                                                         sample_period)
        self._event_gen = EventSpectrumGenerator(np.zeros(rsp.numchans),
                                                 self._sample_per, 
                                                 min_sep=deadtime)

    def set_response(self, rsp):
        """Set/change the detector response.
        
        Args:
            rsp (:class:`~gbm.data.RSP`): A detector response object        
        """
        self._rsp = rsp
        self._spec_gen = VariableSourceSpectrumGenerator(rsp, self._spec_func,
                                                         self._spec_params,
                                                         self._sample_per)

    def set_spectrum(self, spec_func, spec_params):
        """Set/change the spectrum.
        
        Args:
            spec_func (:class:`~gbm.spectra.functions.Function`):
                A photon model function
            spec_params (iterable): The parameters for the function
        """
        self._spec_func = spec_func
        self._spec_params = spec_params
        self._spec_gen = VariableSourceSpectrumGenerator(self._rsp,
                                                         self._spec_func,
                                                         self._spec_params,
                                                         self._sample_per)

    def set_time_profile(self, time_func, time_params):
        """Set/change the time profile.
        
        Args:
            time_func (<function>): A time profile function
            time_params (iterable): Parameters for the time profile function
        """
        self._time_func = time_func
        self._time_params = time_params

    def simulate(self, tstart, tstop):
        """Generate an EventList containing the individual counts from the 
        simulation
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
        
        Returns:
            :class:`~gbm.data.primitives.EventList`:
                The simulated EventList
        """
        # create the time grid
        dur = (tstop - tstart)
        numpts = int(round(dur / self._sample_per))
        time_array = np.linspace(tstart, tstop, numpts)

        # calculate the spectral amplitudes over the grid
        amps = self._time_func(time_array, *self._time_params)

        times = []
        chans = []
        for i in range(numpts):
            # update amplitude and generate the count spectrum
            self._spec_gen.amp = amps[i]
            self._event_gen.spectrum = next(self._spec_gen).counts
            # generate the count arrival times for the time slice spectrum
            events = next(self._event_gen)
            if events is not None:
                times.extend((events[0] + time_array[i]).tolist())
                chans.extend(events[1].tolist())

        # create the eventlist
        eventlist = EventList.from_lists(times, chans,
                                         self._rsp.ebounds['E_MIN'],
                                         self._rsp.ebounds['E_MAX'])
        eventlist.sort('TIME')
        return eventlist

    def to_tte(self, tstart, tstop, trigtime=None, **kwargs):
        """Generate an TTE object containing the individual counts from the 
        simulation
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
            trigtime (float, optional): The trigger time. Default is 0.
            **kwargs: Options to pass to :class:`~gbm.data.TTE`
        
        Returns:
            :class:`~gbm.data.TTE`:
                The simulated TTE object
        """
        if trigtime is None:
            trigtime = 0.0
        eventlist = self.simulate(tstart, tstop)
        tte = TTE.from_data(eventlist, trigtime=trigtime, **kwargs)
        return tte


class TteBackgroundSimulator:
    """Simulate TTE or EventList data given a modeled background spectrum and
    time profile model. The spectrum is fixed throughout the time profile of
    the background, but the amplitude of the background is time-dependent, set
    by the time profile function.

    Parameters:
        bkgd_spectrum (:class:`~gbm.background.BackgroundSpectrum`):
            A modeled background spectrum
        distrib (str): The distribution from which the background is
                       simulated; either 'Poisson' or 'Gaussian'
        time_func (<function>): A time profile function
        time_params (iterable): Parameters for the time profile function
        sample_period (float, optional): The sampling period of the simulator
            in seconds. Default is 0.001. The simulator will produce arrival 
            times consistent with a spectrum over a finite time slice.  This 
            time slice should be short enough to approximate a non-homogeneous 
            Poisson process, but long enough to allow for a tractable 
            computation time.
        deadtime (float, optional): The dead time in seconds for each recorded 
                                    count during which another count cannot be 
                                    recorded. Default is 2e-6 s.
    """
    def __init__(self, bkgd_spectrum, distrib, time_func, time_params,
                 sample_period=0.001, deadtime=2e-6):
        
        warnings.filterwarnings("ignore", category=UserWarning)

        if sample_period <= 0.0:
            raise ValueError('Sample period must be positive')
        if deadtime < 0.0:
            raise ValueError('Deadtime must be non-negative')
        
        self._spec_gen = None
        self._bkgd = bkgd_spectrum
        self._time_func = time_func
        self._time_params = time_params
        self._sample_per = sample_period
        self._deadtime = deadtime
        self._event_gen = EventSpectrumGenerator(np.zeros(self._bkgd.size),
                                                 self._sample_per, 
                                                 min_sep=deadtime)
        self.set_background(bkgd_spectrum, distrib)

    def set_background(self, bkgd_spectrum, distrib):
        """Set/change the spectrum.
        
        Args:
            bkgd_spectrum (:class:`~gbm.background.BackgroundSpectrum`):
                A modeled background spectrum
            distrib (str): The distribution from which the background is
                           simulated; either 'Poisson' or 'Gaussian'
        """
        bkgd = BackgroundSpectrum(bkgd_spectrum.rates,
                                  bkgd_spectrum.rate_uncertainty,
                                  bkgd_spectrum.lo_edges,
                                  bkgd_spectrum.hi_edges,
                                  [self._sample_per] * bkgd_spectrum.size)

        if distrib == 'Poisson':
            self._spec_gen = VariablePoissonBackground(bkgd)
        elif distrib == 'Gaussian':
            self._spec_gen = VariableGaussianBackground(bkgd)
        else:
            raise ValueError("distrib can only be 'Poisson' or 'Gaussian'")

    def simulate(self, tstart, tstop):
        """Generate an EventList containing the individual counts from the 
        background simulation
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
        
        Returns:
            :class:`~gbm.data.primitives.EventList`:
                The simulated EventList
        """
        # create the time grid
        dur = (tstop - tstart)
        numpts = int(round(dur / self._sample_per))
        time_array = np.linspace(tstart, tstop, numpts)

        # calculate the spectral amplitudes over the grid
        amps = self._time_func(time_array, *self._time_params)

        times = []
        chans = []
        for i in range(numpts):
            # update amplitude and generate the count spectrum
            self._spec_gen.amp = amps[i]
            self._event_gen.spectrum = self._whole_counts(
                next(self._spec_gen).counts)
            # generate the count arrival times for the time slice spectrum
            events = next(self._event_gen)
            if events is not None:
                times.extend((events[0] + time_array[i]).tolist())
                chans.extend(events[1].tolist())

        # create the eventlist
        eventlist = EventList.from_lists(times, chans, self._bkgd.lo_edges,
                                         self._bkgd.hi_edges)
        eventlist.sort('TIME')
        return eventlist

    def to_tte(self, tstart, tstop, trigtime=None, **kwargs):
        """Generate an TTE object containing the individual counts from the 
        background simulation
        
        Args:
            tstart (float): The start time of the simulation
            tstop (float): The stop time of the simulation
            trigtime (float, optional): The trigger time. Default is 0.
            **kwargs: Options to pass to :class:`~gbm.data.TTE`
        
        Returns:
            :class:`~gbm.data.TTE`:
                The simulated TTE object
        """
        if trigtime is None:
            trigtime = 0.0
        eventlist = self.simulate(tstart, tstop)
        tte = TTE.from_data(eventlist, trigtime=trigtime, **kwargs)
        return tte

    def _whole_counts(self, counts):
        # because we can end up with fractional counts for the background
        # (the *rate* is what is typically modeled, and so no guarantee that
        #  counts will come out to whole integers)
        u = np.random.random(counts.size)
        whole_counts = counts.astype(int)
        mask = (counts - whole_counts) > u
        whole_counts[mask] += 1
        return whole_counts
