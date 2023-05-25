#! /usr/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin
# Email : shengjunlin@asiaa.sinica.edu.tw
# Description : An example for reading POL2 catlogue created by Starlink
# Date : 2022-03-16
# #########################################################
from astropy.io import fits
from polarization import *

# POL2 catlogue
cat_filename = 'polcat_debias_04arcsec_TYPEas.FIT'

with fits.open(cat_filename, memmap=True) as HDUlist:
    TABdata = HDUlist[1].data
    p = data['P'] / 100.  # polarization fraction
    dp = data['DP'] / 100.  # error of polarization fraction
    pSNR = p/dp  # SNR of polarization fraction

polseg_from_POL2_cat(
        reg_filename='polcat_debias_04arcsec_10per8as_pSNR3_dpLess0.1.reg',
        data        =TABdata,
        mask        =((3. <= pSNR) & (dp <= 0.1)),  # Index slicing
        ten_percent_scale_asec=8.,
        PA_offset   =90.,
        seg_color   ='blue',
        seg_width   =1)

