# scat.py: GBM SCAT file class
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
import astropy.io.fits as fits
from . import headers as hdr
from gbm.time import Met
from .data import DataFile


class Parameter:
    """A fit parameter class
    
    Parameters:
        value (float): The central fit value
        uncert (float or 2-tuple): The 1-sigma uncertainty. If a 2-tuple, then 
                                   is of the form (low, high)
        name (str, optional): The name of the parameter
        units (str, optional): The units of the parameter
        support (2-tuple, optional): The valid support of the parameter
    
    Attributes:
        name (str): The name of the parameter
        support (2-tuple): The valid support of the parameter
        uncertainty (2-tuple): The 1-sigma uncertainty
        units (str): The units of the parameter
        value (float): The central fit value
    """
    def __init__(self, value, uncert, name='', units=None, 
                 support=(-np.inf, np.inf)):
        
        self._value = float(value)
        if isinstance(uncert, (tuple, list)):
            if len(uncert) == 2:
                pass
            elif len(uncert) == 1:
                uncert = (uncert[0], uncert[0])
            else:
                raise ValueError('uncertainty must be a 1- or 2-tuple')
        elif isinstance(uncert, float):
            uncert = (uncert, uncert)
        else:
            raise TypeError('uncertainty must be a float or 1- or 2-tuple')
        self._uncert = uncert
        
        self._units = units
        self._name = name
        self._support = support
    
    @property
    def value(self):
        return self._value
    @property
    def uncertainty(self):
        return self._uncert
    @property
    def name(self):
        return self._name
    @property
    def units(self):
        return self._units
    @property
    def support(self):
        return self._support
    
    def __str__(self):
        value, uncertainty = self._str_format()

        if uncertainty[0] == uncertainty[1]:
            s = '+/- {0}'.format(uncertainty[0])
        else:
            s = '+{0}/-{1}'.format(uncertainty[0], uncertainty[1])
        if self.units is None:
            return '{0}: {1} {2}'.format(self.name, value, s)
        else:
            return '{0}: {1} {2} {3}'.format(self.name, value, s, self.units)

    def _str_format(self):
        if (self.value > 0.005) and (self.uncertainty[0] > 0.005):
            value = '{0:.2f}'.format(self.value)
            uncertainty = tuple(
                ['{0:.2f}'.format(u) for u in self.uncertainty])
        else:
            value = '{0:.2e}'.format(self.value)
            val_coeff, val_exp = value.split('e')
            val_exp = int(val_exp)
            uncertainty = ['{0:.2e}'.format(u) for u in self.uncertainty]
            uncert_coeff = []
            uncert_exp = []
            for uncert in uncertainty:
                uncert_coeff.append(uncert.split('e')[0])
                uncert_exp.append(int(uncert.split('e')[1]))
        return (value, uncertainty)

    def valid_value(self):
        """Check if the parameter value is within the allowed parameter range
        """
        if (self.value >= self.support[0]) and \
                (self.value <= self.support[1]):
            return True
        else:
            return False

    def one_sigma_range(self):
        """Return the 1 sigma range of the parameter fit
        """
        return (self.value - self.uncertainty[0], 
                self.value + self.uncertainty[1])

    def to_fits_value(self):
        """Return as a tuple to be used for a FITS file

        Returns:
            (tuple): 2-value tuple (value, uncertainty) or 3-value tuple 
                     (value, +uncertainty, -uncertainty)   
        """
        return (self.value, *self.uncertainty[::-1])

class PhotonFlux(Parameter):
    """A photon flux. Inherits from :class:`Parameter`.
    
    Parameters:
        value (float): The central flux value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    
    Attributes:
        energy_range (tuple): The enery range (low, high)
        name (str): 'Photon Flux'
        support (2-tuple): (0.0, np.inf)
        uncertainty (2-tuple): The 1-sigma uncertainty
        units (str): 'ph/cm^2/s'
        value (float): The central flux value
    """
    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Photon Flux', units='ph/cm^2/s',
                         support=(0.0, np.inf))
        self._energy_range = energy_range
    
    @property
    def energy_range(self):
        return self._energy_range

class PhotonFluence(Parameter):
    """A photon fluence. Inherits from :class:`Parameter`.
    
    Parameters:
        value (float): The central fluence value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    
    Attributes:
        energy_range (tuple): The enery range (low, high)
        name (str): 'Photon Fluence'
        support (2-tuple): (0.0, np.inf)
        uncertainty (2-tuple): The 1-sigma uncertainty
        units (str): 'ph/cm^2'
        value (float): The central fluence value
    """
    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Photon Fluence', units='ph/cm^2',
                         support=(0.0, np.inf))
        self._energy_range = energy_range
    
    @property
    def energy_range(self):
        return self._energy_range

class EnergyFlux(Parameter):
    """An energy flux. Inherits from :class:`Parameter`.
    
    Parameters:
        value (float): The central flux value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    
    Attributes:
        energy_range (tuple): The enery range (low, high)
        name (str): 'Energy Flux'
        support (2-tuple): (0.0, np.inf)
        uncertainty (2-tuple): The 1-sigma uncertainty
        units (str): 'erg/cm^2/s'
        value (float): The central flux value
    """
    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Energy Flux', units='erg/cm^2/s',
                         support=(0.0, np.inf))
        self._energy_range = energy_range
    
    @property
    def energy_range(self):
        return self._energy_range

class EnergyFluence(Parameter):
    """An energy fluence. Inherits from :class:`Parameter`.
    
    Parameters:
        value (float): The central fluence value
        uncert (float or 2-tuple): The 1-sigma uncertainty
        energy_range (tuple): A 2-tuple (low, high) for the energy range
    
    Attributes:
        energy_range (tuple): The enery range (low, high)
        name (str): 'Energy Fluence'
        support (2-tuple): (0.0, np.inf)
        uncertainty (2-tuple): The 1-sigma uncertainty
        units (str): 'erg/cm^2'
        value (float): The central fluence value
    """
    def __init__(self, value, uncert, energy_range):
        super().__init__(value, uncert, name='Energy Fluence', units='erg/cm^2',
                         support=(0.0, np.inf))
        self._energy_range = energy_range
    
    @property
    def energy_range(self):
        return self._energy_range


class ModelFit:
    """A container for the info from a model fit
    
    Parameters:
        name (str): The name of the model
        time_range (float, float): The time range of the model fit, (low, high)
        parameters (list, optional): A list of model parameters
        photon_flux (:class:`PhotonFlux`, optional): The photon flux
        energy_flux (:class:`EnergyFlux`, optional): The energy flux
        photon_fluence (:class:`PhotonFluence`, optional): The photon fluence
        energy_fluence (:class:`EnergyFluence`, optional): The energy fluence
        flux_energy_range (tuple, optional): The energy range of the flux
                                             and fluence, (low, high)
        stat_name (str, optional): The name of the fit statistic
        stat_value (float, optional): The fit statistic value
        dof (int, optional): The degrees-of-freedom of the fit
        covariance (np.array, optional): The covariance matrix of the fit

    Attributes:
        covariance (np.array): The covariance matrix of the fit
        dof (int): The degrees-of-freedom of the fit
        energy_fluence (:class:`EnergyFluence`): The energy fluence
        energy_flux (:class:`EnergyFlux``): The energy flux
        flux_energy_range (tuple): The energy range of the flux and fluence, 
                                   (low, high)
        name (str): The name of the model
        parameters (list): A list of model parameters
        photon_fluence (:class:`PhotonFluence`): The photon fluence
        photon_flux (:class:`PhotonFlux`): The photon flux
        stat_name (str): The name of the fit statistic
        stat_value (float): The fit statistic value
        time_range (float, float): The time range of the model fit, (low, high)
    """
    def __init__(self, name, time_range, **kwargs):
        self._name = str(name)
        if not isinstance(time_range, (list, tuple)):
            raise ValueError('time_range must be a 2-tuple')
        else:
            if len(time_range) != 2:
                raise ValueError('time_range must be a 2-tuple')
        self._time_range = time_range
        
        self._parameters = []
        self._photon_flux = None
        self._energy_flux = None
        self._photon_fluence = None
        self._energy_fluence = None
        self._flux_energy_range = None
        self._stat_name = None
        self._stat_value = None
        self._dof = None
        self._covariance = None

        self._init_by_dict(kwargs)

    def __str__(self):
        param_str = '\n   '.join([str(param) for param in self.parameters])
        return '{0}\n   {1}'.format(self.name, param_str)

    @property
    def name(self):
        return self._name
    @property
    def time_range(self):
        return self._time_range
    
    @property
    def parameters(self):
        return self._parameters
    @parameters.setter
    def parameters(self, val):
        if not isinstance(val, (list, tuple)):
            raise TypeError('parameters must be a list of parameters')
        for p in val:
            if not isinstance(p, Parameter):
                raise TypeError('parameters must be of Parameter type')
        self._parameters = val
    
    @property
    def photon_flux(self):
        return self._photon_flux
    @photon_flux.setter
    def photon_flux(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_flux must be of Parameter type')
        self._photon_flux = val

    @property
    def energy_flux(self):
        return self._energy_flux
    @energy_flux.setter
    def energy_flux(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_flux must be of Parameter type')
        self._energy_flux = val

    @property
    def photon_fluence(self):
        return self._photon_fluence
    @photon_fluence.setter
    def photon_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_fluence must be of Parameter type')
        self._photon_fluence = val

    @property
    def energy_fluence(self):
        return self._energy_fluence
    @energy_fluence.setter
    def energy_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_fluence must be of Parameter type')
        self._energy_fluence = val
    
    @property
    def flux_energy_range(self):
        return self._flux_energy_range
    @flux_energy_range.setter
    def flux_energy_range(self, val):
        if not isinstance(val, (list, tuple)):
            raise ValueError('flux_energy_range must be a 2-tuple')
        else:
            if len(val) != 2:
                raise ValueError('flux_energy_range must be a 2-tuple')
        self._flux_energy_range = val
    
    @property
    def stat_name(self):
        return self._stat_name
    @stat_name.setter
    def stat_name(self, val):
        self._stat_name = str(val)
    
    @property
    def stat_value(self):
        return self._stat_value
    @stat_value.setter
    def stat_value(self, val):
        try:
            float_val = float(val)
        except:
            raise TypeError('stat_value must be a float')
        self._stat_value = float_val
    
    @property
    def dof(self):
        return self._dof
    @dof.setter
    def dof(self, val):
        try:
            int_val = int(val)
        except:
            raise TypeError('dof must be an integer')
        self._dof = int_val
    
    @property
    def covariance(self):
        return self._covariance
    @covariance.setter
    def covariance(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('covariance must an array')
        if len(val.shape) != 2:
            raise ValueError('covariance must be a n x n array')
        if val.shape[0] != val.shape[1]:
            raise ValueError('covariance must be a n x n array')
        
        self._covariance = val
    
    def parameter_list(self):
        """Return the list of parameter names
        
        Returns:
            (list): The parameter names
        """
        return [param.name for param in self.parameters]

    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
            (astropy.io.fits.BinTableHDU): The FITS table
        """
        numparams = len(self.parameters)
        cols = []
        cols.append(
            fits.Column(name='TIMEBIN', format='2D', array=[self.time_range]))
        i = 0
        for param in self.parameters:
            col = fits.Column(name='PARAM{0}'.format(i), format='3E',
                              array=[param.to_fits_value()])
            cols.append(col)
            i += 1

        cols.append(fits.Column(name='PHTFLUX', format='3E',
                                array=[self.photon_flux.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLNC', format='3E',
                                array=[self.photon_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLUX', format='3E',
                                array=[self.energy_flux.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNC', format='3E',
                                array=[self.energy_fluence.to_fits_value()]))
        cols.append(fits.Column(name='REDCHSQ', format='2E',
                                array=[self.stat_value / self.dof]))
        cols.append(fits.Column(name='CHSQDOF', format='1I', array=[self.dof]))
        cols.append(fits.Column(name='COVARMAT',
                                format='{0}E'.format(numparams * numparams),
                                dim='({0},{0})'.format(numparams),
                                array=[self.covariance]))

        hdu = fits.BinTableHDU.from_columns(cols, name='FIT PARAMS')
        return hdu
    
    def _init_by_dict(self, values):
        for key, val in values.items():
            try:
                p = getattr(self, '_'+key)
                if isinstance(p, property):
                    getattr(self, key).__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute".format(key))

class GbmModelFit(ModelFit):
    """A container for the info from a model fit, with values used in the 
    GBM SCAT files. Inherits from :class:`ModelFit`.

    Attributes:
        covariance (np.array): The covariance matrix of the fit
        dof (int): The degrees-of-freedom of the fit
        duration_fluence (:class:`EnergyFluence`): The energy fluence over the
                                                   duration energy range, 
                                                   nominally 50-300 keV
        energy_fluence (:class:`EnergyFluence`): The energy fluence, nominally 
                                                 over 10-1000 keV
        energy_fluence_50_300 (:class:`EnergyFluence`): The energy fluence over
                                                        50-300 keV
        energy_flux (:class:`EnergyFlux`): The energy flux, nominally over 
                                           10-1000 keV
        flux_energy_range (tuple): The energy range of the flux and fluence, 
                                   (low, high)
        name (str): The name of the model
        parameters (list): A list of model parameters
        photon_fluence (:class:`PhotonFluence`): The photon fluence, nominally 
                                                 over 10-1000 keV
        photon_flux (:class:`PhotonFlux`): The photon flux, nominally over 
                                           10-1000 keV
        photon_flux_50_300 (:class:`PhotonFlux`): The photon flux over 50-300 keV
        stat_name (str): The name of the fit statistic
        stat_value (float): The fit statistic value
        time_range (float, float): The time range of the model fit, (low, high)
    """
    def __init__(self, name, time_range, **kwargs):
        
        self._photon_flux_50_300 = None
        self._energy_fluence_50_300 = None
        self._duration_fluence = None
    
        super().__init__(name, time_range, **kwargs)
    
    @property
    def photon_flux_50_300(self):
        return self._photon_flux_50_300
    @photon_flux_50_300.setter
    def photon_flux_50_300(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_flux_50_300 must be of Parameter type')
        self._photon_flux_50_300 = val

    @property
    def energy_fluence_50_300(self):
        return self._energy_fluence_50_300
    @energy_fluence_50_300.setter
    def energy_fluence_50_300(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_fluence_50_300 must be of Parameter type')
        self._energy_fluence_50_300 = val

    @property
    def duration_fluence(self):
        return self._duration_fluence
    @duration_fluence.setter
    def duration_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('duration_fluence must be of Parameter type')
        self._duration_fluence = val

    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
            (astropy.io.fits.BinTableHDU): The FITS table
        """
        numparams = len(self.parameters)
        cols = []
        cols.append(
            fits.Column(name='TIMEBIN', format='2D', array=[self.time_range]))
        i = 0
        for param in self.parameters:
            col = fits.Column(name='PARAM{0}'.format(i), format='3E',
                              array=[param.to_fits_value()])
            cols.append(col)
            i += 1

        cols.append(fits.Column(name='PHTFLUX', format='3E',
                                array=[self.photon_flux.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLNC', format='3E',
                                array=[self.photon_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLUX', format='3E',
                                array=[self.energy_flux.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNC', format='3E',
                                array=[self.energy_fluence.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLUXB', format='3E',
                                array=[
                                    self.photon_flux_50_300.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNCB', format='3E',
                                array=[
                                    self.energy_fluence_50_300.to_fits_value()]))
        cols.append(fits.Column(name='DURFLNC', format='3E',
                                array=[self.duration_fluence.to_fits_value()]))
        cols.append(fits.Column(name='REDCHSQ', format='2E',
                                array=[[self.stat_value / self.dof] * 2]))
        cols.append(fits.Column(name='CHSQDOF', format='1I', array=[self.dof]))
        cols.append(fits.Column(name='COVARMAT',
                                format='{0}E'.format(numparams * numparams),
                                dim='({0},{0})'.format(numparams),
                                array=[self.covariance]))

        hdu = fits.BinTableHDU.from_columns(cols, name='FIT PARAMS')
        return hdu

    @classmethod
    def from_fits_row(cls, fits_row, model_name, param_names=None,
                      flux_range=(10.0, 1000.0), dur_range=(50.0, 300.0)):
        """Read a FITS row and return a :class:`GbmModelFit` object
        
        Returns:
            (:class:`GbmModelFit`)
        """
        time_range = tuple(fits_row['TIMEBIN'])
        nparams = sum([1 for name in fits_row.array.dtype.names \
                       if 'PARAM' in name])
        if param_names is None:
            param_names = ['']*nparams
        
        params = []
        for i in range(nparams):
            param = fits_row['PARAM' + str(i)]
            params.append(Parameter(param[0], tuple(param[1:]), 
                          name=param_names[i]))
        
        pflux = PhotonFlux(fits_row['PHTFLUX'][0], 
                           tuple(fits_row['PHTFLUX'][1:]), flux_range)
        pflnc = PhotonFluence(fits_row['PHTFLNC'][0],
                              tuple(fits_row['PHTFLNC'][1:]), flux_range)
        eflux = EnergyFlux(fits_row['NRGFLUX'][0],
                           tuple(fits_row['NRGFLUX'][1:]), flux_range)
        eflnc = EnergyFluence(fits_row['NRGFLNC'][0],
                              tuple(fits_row['NRGFLNC'][1:]), flux_range)
        pfluxb = PhotonFlux(fits_row['PHTFLUXB'][0], 
                            tuple(fits_row['PHTFLUXB'][1:]), (50.0, 300.0))
        eflncb = EnergyFluence(fits_row['NRGFLNCB'][0], 
                               tuple(fits_row['NRGFLNCB'][1:]), (50.0, 300.0))
        durflnc = PhotonFluence(fits_row['DURFLNC'][0],
                                tuple(fits_row['DURFLNC'][1:]), dur_range)
        dof = fits_row['CHSQDOF']
        # scat provides the [fit stat, chisq], while bcat is only the fit stat
        try:
            stat_val = fits_row['REDCHSQ'][0]*dof
        except:
            stat_val = fits_row['REDCHSQ']*dof
        covar = fits_row['COVARMAT']
                             
        obj = cls(model_name, time_range, parameters=params, photon_flux=pflux,
                  photon_fluence=pflnc, energy_flux=eflux, energy_fluence=eflnc,
                  flux_energy_range=flux_range, stat_value=stat_val, dof=dof,
                  covariance=covar, photon_flux_50_300=pfluxb, 
                  energy_fluence_50_300=eflncb, duration_fluence=durflnc)
        return obj
        

class DetectorData():
    """A container for detector info used in a fit
    
    Parameters:
        instrument (str): The name of the instrument
        detector (str): The name of the detector
        datatype (str): The name of the datatype
        filename (str): The filename of the data file
        numchans (int): Number of energy channels used
        active (bool, optional): True if the detector is used in the fit
        response (str, optional): The filename of the detector response
        time_range (tuple, optional): The time range of the data used
        energy_range (tuple, optional): The energy range of the data used
        channel_range (tuple, optional): The energy channel range of the data
        energy_edges (np.array, optional): The edges of the energy channels
        photon_counts (np.array, optional): The deconvolved photon counts for 
                                            the detector
        photon_model (np.array, optional): The photon model for the detector
        photon_errors (np.array, optional): The deconvolved photon count errors 
                                            for the detector

    Attributes:
        active (bool, optional): True if the detector is used in the fit
        channel_range = (int, int): The energy channel range of the data
        datatype (str): The name of the datatype
        detector (str): The name of the detector
        energy_edges (np.array): The edges of the energy channels
        energy_range (float, float): The energy range of the data used
        filename (str): The filename of the data file
        instrument (str): The name of the instrument
        numchans (int): Number of energy channels used
        photon_counts (np.array): The deconvolved photon counts for the detector
        photon_errors (np.array): The deconvolved photon count errors for the 
                                  detector
        photon_model (np.array): The photon model for the detector
        response (str): The filename of the detector response
        time_range (float, float): The time range of the data used
    """
    def __init__(self, instrument, detector, datatype, filename, numchans,
                 **kwargs):
        self._instrument = instrument
        self._detector = detector
        self._datatype = datatype
        self._filename = filename
        self._numchans = int(numchans)
        
        self._active = True
        self._response = ''
        self._time_range = (None, None)
        self._energy_range = (None, None)
        self._channel_range = (None, None)
        self._channel_mask = None
        self._energy_edges = None
        self._photon_counts = None
        self._photon_model = None
        self._photon_errors = None
        
        self._init_by_dict(kwargs)
            
    # read-only
    @property
    def instrument(self):
        return self._instrument
    @property
    def detector(self):
        return self._detector
    @property
    def datatype(self):
        return self._datatype
    @property
    def filename(self):
        return self._filename
    @property
    def numchans(self):
        return self._numchans

    @property
    def active(self):
        return self._active
    @active.setter
    def active(self, val):
        try:
            bool_val = bool(val)
        except:
            raise TypeError('active must be Boolean')
        self._active = bool_val

    @property
    def response(self):
        return self._response
    @response.setter
    def response(self, val):
        if not isinstance(val, str):
            raise TypeError('response filename must be a string')
        self._response = val

    @property
    def time_range(self):
        return self._time_range
    @time_range.setter
    def time_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('time_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('time_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('time_range must be of form (low, high)')
        else:
            pass
        
        self._time_range = val

    @property
    def energy_range(self):
        return self._energy_range
    @energy_range.setter
    def energy_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('energy_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('energy_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('energy_range must be of form (low, high)')
        else:
            pass
        
        self._energy_range = val
        
    @property
    def channel_range(self):
        return self._channel_range
    @channel_range.setter
    def channel_range(self, val):
        if not isinstance(val, (tuple, list)):
            raise TypeError('channel_range must be a 2-tuple')
        elif len(val) != 2:
            raise ValueError('channel_range must be a 2-tuple')
        elif val[0] > val[1]:
            raise ValueError('channel_range must be of form (low, high)')
        else:
            pass
        
        self._channel_range = val

    @property
    def channel_mask(self):
        return self._channel_mask
    @channel_mask.setter
    def channel_mask(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('channel_mask must be an array')
        self._channel_mask = val.astype(bool)

    @property
    def energy_edges(self):
        return self._energy_edges
    @energy_edges.setter
    def energy_edges(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('energy_edges must be an array')
        self._energy_edges = val

    @property
    def photon_counts(self):
        return self._photon_counts
    @photon_counts.setter
    def photon_counts(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_counts must be an array')
        self._photon_counts = val

    @property
    def photon_model(self):
        return self._photon_model
    @photon_model.setter
    def photon_model(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_model must be an array')
        self._photon_model = val

    @property
    def photon_errors(self):
        return self._photon_errors
    @photon_errors.setter
    def photon_errors(self, val):
        if not isinstance(val, np.ndarray):
            raise TypeError('photon_errors must be an array')
        self._photon_errors = val

    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
            (astropy.io.fits.BinTableHDU): The FITS row
        """
        numchans = len(self.energy_edges)
        e_dim = str(numchans) + 'E'
        p_dim = str(numchans - 1) + 'E'
        p_unit = 'Photon cm^-2 s^-1 keV^-1'
        fit_int = '{0}: {1} s, '.format(self.time_range[0], self.time_range[1])
        fit_int += '{0}: {1} keV, '.format(self.energy_range[0],
                                           self.energy_range[1])
        fit_int += 'channels {0}: {1}'.format(self.channel_range[0],
                                              self.channel_range[1])

        col1 = fits.Column(name='INSTRUME', format='20A',
                           array=[self.instrument])
        col2 = fits.Column(name='DETNAM', format='20A', array=[self.detector])
        col3 = fits.Column(name='DATATYPE', format='20A',
                           array=[self.datatype])
        col4 = fits.Column(name='DETSTAT', format='20A',
                           array=['INCLUDED' if self.active else 'OMITTED'])
        col5 = fits.Column(name='DATAFILE', format='60A',
                           array=[self.filename])
        col6 = fits.Column(name='RSPFILE', format='60A', array=[self.response])
        col7 = fits.Column(name='FIT_INT', format='60A', array=[fit_int])
        col8 = fits.Column(name='CHANNUM', format='1J', array=[numchans - 1])
        col9 = fits.Column(name='FITCHAN', format='{}J'.format(numchans),
                           array=[self.channel_mask])
        col10 = fits.Column(name='E_EDGES', format=e_dim, unit='keV',
                            array=[self.energy_edges])
        col11 = fits.Column(name='PHTCNTS', format=p_dim, unit=p_unit,
                            array=[self.photon_counts])
        col12 = fits.Column(name='PHTMODL', format=p_dim, unit=p_unit,
                            array=[self.photon_model])
        col13 = fits.Column(name='PHTERRS', format=p_dim, unit=p_unit,
                            array=[self.photon_errors])
        hdu = fits.BinTableHDU.from_columns(
            [col1, col2, col3, col4, col5, col6,
             col7, col8, col9, col10, col11,
             col12, col13], name='DETECTOR DATA')
        return hdu

    @classmethod
    def from_fits_row(cls, fits_row):
        """Read a FITS row and return a DetectorData object
        
        Returns:
            (:class:`DetectorData`)
        """
        instrument = fits_row['INSTRUME']
        det = fits_row['DETNAM']
        datatype = fits_row['DATATYPE']
        if fits_row['DETSTAT'] == 'INCLUDED':
            active = True
        else:
            active = False
        datafile = fits_row['DATAFILE']
        rspfile = fits_row['RSPFILE']
        fit_ints = fits_row['FIT_INT'].split(' ')
        time_range = (float(fit_ints[0][:-1]), float(fit_ints[1]))
        energy_range = (float(fit_ints[3][:-1]), float(fit_ints[4]))
        channel_range = (int(fit_ints[8][:-1]), int(fit_ints[9]))
        numchans = fits_row['CHANNUM']
        channel_mask = np.zeros(numchans, dtype=bool)
        channel_mask[fits_row['FITCHAN'][0]:fits_row['FITCHAN'][1]+1] = True
        energy_edges = fits_row['E_EDGES']
        photon_counts = fits_row['PHTCNTS']
        photon_model = fits_row['PHTMODL']
        photon_errors = fits_row['PHTERRS']             
        obj = cls(instrument, det, datatype, datafile, numchans, active=active,
                  response=rspfile, time_range=time_range, 
                  energy_range=energy_range, channel_range=channel_range, 
                  channel_mask=channel_mask, energy_edges=energy_edges, 
                  photon_counts=photon_counts, photon_model=photon_model, 
                  photon_errors=photon_errors)
        return obj

    def _init_by_dict(self, values):
        for key, val in values.items():
            try:
                p = getattr(self, '_'+key)
                if isinstance(p, property):
                    getattr(self, key).__set__(self, val)
                else:
                    self.__setattr__(key, val)
            except AttributeError:
                raise ValueError("{} is not a valid attribute".format(key))
    

class Scat(DataFile):
    """A container class for the spectral fit data in an SCAT file
    
    Attributes:
        detectors (list): The :class:`DetectorData` objects used in the analysis
        headers (dict): The SCAT file headers
        model_fits (list): The :class:`GbmModelFit` objects, one for each model 
                           fit
        num_detectors (int): The number of detectors in the SCAT file
        num_fits (int): The number of model fits
    """
    def __init__(self):
        self._detectors = []
        self._model_fits = []
        self._headers = {}
    
    @property
    def detectors(self):
        return self._detectors
    @property
    def model_fits(self):
        return self._model_fits
    @property
    def headers(self):
        return self._headers
    @property
    def num_detectors(self):
        return len(self.detectors)
    @property
    def num_fits(self):
        return len(self.model_fits)
    
    def add_detector_data(self, detector_data):
        """Add a new detector to the Scat
        
        Args:
            detector_data (:class:`DetectorData`): The detector data
        """
        if not isinstance(detector_data, DetectorData):
            raise TypeError("Can only add DetectorData objects")
        self._detectors.append(detector_data)

    def add_model_fit(self, model_fit):
        """Add a new model fit to the Scat
        
        Args:
            model_fit (:class:`GbmModelFit`): The model fit data
        """
        if not isinstance(model_fit, GbmModelFit):
            raise TypeError("Can only add GbmModelFit objects")
        self._model_fits.append(model_fit)

    @classmethod
    def open(cls, filename):
        """Open a SCAT FITS file and create a Scat object
        
        Args:
            filename (str): The file to open
        
        Returns:
            (:class:`Scat`)
        """
        obj = cls()
        obj._file_properties(filename)
        # open FITS file
        with fits.open(filename) as hdulist:
            for hdu in hdulist:
                obj._headers.update({hdu.name: hdu.header})

            # read the detector data HDU
            det_data = hdulist['DETECTOR DATA'].data
            for row in det_data:
                obj.add_detector_data(DetectorData.from_fits_row(row))

            # read the fit params HDU
            fit_hdu = hdulist['FIT PARAMS']
            obj._from_fitparam_hdu(fit_hdu)
            
        return obj

    def write(self, directory, filename=None):
        """Write a Scat object to a FITS file
        
        Args:
            directory (str): The directory where the file is to be written
            filename (str, optional): The filename.  If not set, a default 
                                      filename will be used.
        """
        raise NotImplementedError
        
        if (self.filename is None) and (filename is None):
            raise NameError('Filename not set')
        if filename is None:
            filename = self.filename_obj.basename()
        self.set_filename(filename, directory=directory)

        # initialize the FITS file
        hdulist = fits.HDUList()
        prihdr = self._primary_header()
        primary_hdu = fits.PrimaryHDU(header=prihdr)
        primary_hdu.add_checksum()
        hdulist.append(primary_hdu)

        # construct the detector data extension
        det_hdu = None
        for detector in self._detectors:
            if det_hdu is None:
                det_hdu = detector.to_fits_row()
            else:
                det_hdu.data = np.concatenate(
                    (det_hdu.data, detector.to_fits_row().data))
        det_hdu.header = self._update_detector_hdr(det_hdu.header)
        det_hdu.add_checksum()
        hdulist.append(det_hdu)

        # construct the fit extension
        fit_hdu = None
        for fit in self._model_fits:
            if fit_hdu is None:
                fit_hdu = fit.to_fits_row()
            else:
                fit_hdu.data = np.concatenate(
                    (fit_hdu.data, fit.to_fits_row().data))
        fit_hdu.header = self._update_fitparam_hdr(fit_hdu.header)
        fit_hdu.add_checksum()
        hdulist.append(fit_hdu)

        # write out the file
        filename = directory + filename
        hdulist.writeto(self.filename, checksum=True, clobber=True)

    def _from_fitparam_hdu(self, fit_hdu):

        fit_data = fit_hdu.data
        fit_hdr = fit_hdu.header
        nparams = fit_hdr['N_PARAM']
        
        # find unique models in the event there are multiple components
        models = [fit_hdr.comments['TTYPE' + str(2 + i)].split(':')[0]
                  for i in range(nparams)]
        model = '+'.join(list(set(models)))

        # the parameter names
        param_names = [
            fit_hdr.comments['TTYPE' + str(2 + i)].split(':')[1].strip()
            for i in range(nparams)]
        
        # populate each model fit
        for row in fit_data:
            modelfit = GbmModelFit.from_fits_row(row, model,
                                                 param_names=param_names)
            modelfit.stat_name = fit_hdr['STATISTC']
            self.add_model_fit(modelfit)

    def _update_detector_hdr(self, det_hdr):
        det_hdr.comments['TTYPE1'] = 'Instrument name for this detector'
        det_hdr.comments[
            'TTYPE2'] = 'Detector number; if one of several available'
        det_hdr.comments['TTYPE3'] = 'Data type used for this analysis'
        det_hdr.comments['TTYPE4'] = 'Was this detector INCLUDED or OMITTED'
        det_hdr.comments['TTYPE5'] = 'Data file name for this dataset'
        det_hdr.comments['TTYPE6'] = 'Response file name for this dataset'
        det_hdr.comments['TTYPE7'] = 'Fit intervals'
        det_hdr.comments[
            'TTYPE8'] = 'Total number of energy channels for this detector'
        det_hdr.comments[
            'TTYPE9'] = 'Channels selected in fitting this detector'
        det_hdr.comments['TTYPE10'] = 'Energy edges for each selected detector'
        det_hdr.comments['TTYPE11'] = 'Array of photon counts data'
        det_hdr.comments['TTYPE12'] = 'Array of photon model data'
        det_hdr.comments['TTYPE13'] = 'Array of errors in photon counts data'
        det_hdr['NUMFITS'] = (
            len(self._model_fits), 'Number of spectral fits in the data')
        prihdu = self._primary_header()
        keys = ['ORIGIN', 'TELESCOP', 'INSTRUME', 'OBSERVER', 'MJDREFI',
                'MJDREFF',
                'TIMESYS', 'TIMEUNIT', 'DATE-OBS', 'DATE-END', 'TSTART',
                'TSTOP',
                'TRIGTIME']
        for key in keys:
            det_hdr[key] = (prihdu[key], prihdu.comments[key])
        return det_hdr

    def _update_fitparam_hdr(self, fit_hdr):
        e_range = self._model_fits[0].flux_energy_range
        model_name = self._model_fits[0].name
        param_names = self._model_fits[0].parameter_list()
        numparams = len(param_names)
        statistic = self._model_fits[0].stat_name

        g_range = '({0}-{1} keV)'.format(e_range[0], e_range[1])
        b_range = '(50-300 keV)'

        fit_hdr.comments['TTYPE1'] = 'Start and stop times relative to trigger'
        for i in range(numparams):
            colname = 'TTYPE{0}'.format(i + 2)
            fit_hdr.comments[colname] = '{0}: {1}'.format(model_name,
                                                          param_names[i])
        ttypes = ['TTYPE' + str(numparams + 2 + i) for i in range(10)]
        fit_hdr.comments[
            ttypes[0]] = 'Photon Flux (ph/s-cm^2) std energy ' + g_range
        fit_hdr.comments[
            ttypes[1]] = 'Photon Fluence (ph/cm^2) std energy ' + g_range
        fit_hdr.comments[
            ttypes[2]] = 'Energy Flux (erg/s-cm^2) std energy ' + g_range
        fit_hdr.comments[
            ttypes[3]] = 'Energy Fluence (erg/cm^2) std energy ' + g_range
        fit_hdr.comments[
            ttypes[4]] = 'Reduced Chi^2 (1) and fitting statistic (2)'
        fit_hdr.comments[ttypes[5]] = 'Degrees of Freedom'
        fit_hdr.comments[
            ttypes[6]] = 'Photon Flux (ph/s-cm^2) BATSE energy ' + b_range
        fit_hdr.comments[
            ttypes[7]] = 'Photon Fluence (ph/cm^2) for durations (user)'
        fit_hdr.comments[
            ttypes[8]] = 'Energy Fluence (erg/cm^2) BATSE energy ' + b_range
        fit_hdr.comments[
            ttypes[9]] = 'Covariance matrix for the fir (N_PARAM^2)'
        fit_hdr['N_PARAM'] = (
            numparams, 'Total number of fit parameters (PARAMn)')
        fit_hdr['FLU_LOW'] = (
            e_range[0], 'Lower limit of flux/fluence integration (keV)')
        fit_hdr['FLU_HIGH'] = (
            e_range[1], 'Upeer limit of flux/fluence integration (keV)')
        fit_hdr['STATISTC'] = (
            statistic, 'Indicates merit function used for fitting')
        fit_hdr['NUMFITS'] = (
            len(self._model_fits), 'Number of spectral fits in the data')

        prihdu = self._primary_header()
        keys = ['ORIGIN', 'TELESCOP', 'INSTRUME', 'OBSERVER', 'MJDREFI',
                'MJDREFF',
                'TIMESYS', 'TIMEUNIT', 'DATE-OBS', 'DATE-END', 'TSTART',
                'TSTOP',
                'TRIGTIME']
        for key in keys:
            fit_hdr[key] = (prihdu[key], prihdu.comments[key])
        return fit_hdr

    def _primary_header(self):
        # create standard GBM primary header
        filetype = 'SPECTRAL FITS'
        prihdr = hdr.primary(filetype=filetype, tstart=self.tstart,
                             filename=self.filename_obj.basename(),
                             tstop=self.tstop, trigtime=self.trigtime)
        # remove the keywords we don't need
        del prihdr['DETNAM'], prihdr['OBJECT'], prihdr['RA_OBJ'], \
            prihdr['DEC_OBJ'], prihdr['ERR_RAD'], prihdr['INFILE01']

        return prihdr
