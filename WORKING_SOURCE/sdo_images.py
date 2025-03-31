#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 13:18:45 2025

@author: coletamburri
"""

import aiapy
import astropy
import astropy.time
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sunpy
import sunpy.map
#from aiapy.calibrate.util import get_correction_table
#from aiapy.psf import deconvolve, psf
#from aiapy.response import Channel
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, LogStretch, time_support
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time

# Increases the figure size in this notebook.
mpl.rcParams["savefig.dpi"] = 400
mpl.rcParams["figure.dpi"] = 400

t_start = parse_time("2024-08-08T19:00:00")
search_results = Fido.search(
    a.Time(t_start, t_start + 12 * u.s),
    a.Instrument.aia,
    a.Wavelength(171 * u.angstrom) | a.Wavelength(1600 * u.angstrom) | a.Wavelength(304 * u.angstrom) | a.Wavelength(131 * u.angstrom),
)

files = Fido.fetch(search_results, max_conn=1)

m_171 = sunpy.map.Map(files[0])
m_171.peek(vmin=0)