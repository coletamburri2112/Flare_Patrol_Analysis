# profiles.py: Functions for lightcurve and background time profiles
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
from .generators import *


# pulse shapes
def tophat(x, amp, tstart, tstop):
    """A tophat (rectangular) pulse function.
    
    Args:
        x (np.array): Array of times
        amp (float): The tophat amplitude
        tstart (float): The start time of the tophat
        tstop (float): The end time of the tophat
    
    Returns:
        np.array: The tophat evaluated at ``x`` times
    """
    mask = (x >= tstart) & (x <= tstop)
    fxn = np.zeros_like(x)
    fxn[mask] = amp
    return fxn


def norris(x, amp, tstart, t_rise, t_decay):
    r"""A Norris pulse-shape function:

    :math:`I(t) = A \lambda e^{-\tau_1/t - t/\tau_2} \text{ for } t > 0;\\ 
    \text{ where } \lambda = e^{2\sqrt(\tau_1/\tau_2)};`
    
    and where
    
    * :math:`A` is the pulse amplitude
    * :math:`\tau_1` is the rise time
    * :math:`\tau_2` is the decay time
    
    References:
        `Norris, J. P., et al. 2005 ApJ 627 324
        <https://iopscience.iop.org/article/10.1086/430294>`_
    
    Args:
        x (np.array): Array of times
        amp (float): The amplitude of the pulse
        tstart (float): The start time of the pulse
        t_rise (float): The rise timescal of the pulse
        t_decay (flaot): The decay timescale of the pulse
    
    Returns:
        np.array: The Norris pulse shape evaluated at ``x`` times
    """
    x = np.asarray(x)
    fxn = np.zeros_like(x)
    mask = (x > tstart)
    lam = amp * np.exp(2.0 * np.sqrt(t_rise / t_decay))
    fxn[mask] = lam * np.exp(
        -t_rise / (x[mask] - tstart) - (x[mask] - tstart) / t_decay)
    return fxn


# ------------------------------------------------------------------------------

# background profiles
def constant(x, amp):
    """A constant background function.
    
    Args:
        x (np.array): Array of times
        amp (float): The background amplitude
    
    Returns:
        np.array: The background evaluated at ``x`` times
    """
    fxn = np.empty(x.size)
    fxn.fill(amp)
    return fxn


def linear(x, c0, c1):
    """A linear background function.
    
    Args:
        x (np.array): Array of times
        c0 (float): The constant coefficient
        c1 (float): The linear coefficient
    
    Returns:
        np.array: The background evaluated at ``x`` times
    """
    fxn = c0 + c1 * x
    return fxn


def quadratic(x, c0, c1, c2):
    """A quadratic background function.
    
    Args:
        x (np.array): Array of times
        c0 (float): The constant coefficient
        c1 (float): The linear coefficient
        c2 (float): The quadratic coefficient
    
    Returns:
        np.array: The background evaluated at ``x`` times
    """
    fxn = linear(x, c0, c1) + c2 * x ** 2
    return fxn
