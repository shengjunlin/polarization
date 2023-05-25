#! /usr/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin
# Email : shengjunlin@asiaa.sinica.edu.tw
# Description : An example for reading images created by CASA
# Date : 2022-03-16
# #########################################################
from polarization import *

def mask(i, j, RA_deg, Dec_deg, I, PI, PA_deg):
    """An example of the masking function.

    Please don't modify the sequence of the input parameters of the masking function.
    The above parameters are required because of the design of the function `polseg_from_2darray'.

    Parameters
    ----------
    i : int
        Pixel value for x.
    j : int
        Pixel value for y.
    RA_deg : float
        RA value in degree.
    Dec_deg : float
        Dec value in degree.
    I : float
        Total intensity value in the given unit.
    PI : float
        Polarized intensity value in the given unit.
    PA_deg : float
        Polarization angle value in degree.

    This is an example for filtering pixels with 10- and 3-sigma detection of
    total intensities (the input I) and polarized intensities (the input PI), respectively.
    """
    # User can modify the following criteria.
    rms_I = 0.
    rms_PI = 0.
    p = PI/I  # polarization fraction

    if (I >= 10*rms_I) and (PI >= 3*rms_PI):
        return True   # Create a ds9 segment region
    else:
        return False  # Don't create a ds9 region

# Created by the casa command `immath'.
I_map_filename      = '.fits'
PI_map_filename     = '.fits'
PA_map_deg_filename = '.fits'

polseg_from_images(reg_filename,
        I_map=I_map_filename,
        PI_map=PI_map_filename,
        PA_map_deg=PA_map_deg_filename,
        i_hdu=0,
        i_chan=0,  # Always 0 for the continuum polarization data.
        sampling_interval_px=3,
        ten_percent_scale_asec=0.05,
        uniform_scale_asec=None,
        uniform_PA_deg=None,
        PA_offset=0.,
        I_clip=0.,  # In the given unit. This can be replaced by `mask_func'.
        PI_clip=0., # In the given unit. This can be replaced by `mask_func'.
        mask_func=None,  # `mask_func=mask' if the user wants to use the above mask function.
        seg_color='black',
        seg_width=1,
        frame='icrs',
        hist_plot=True)  # Show a histrogram in the interactive mode.
