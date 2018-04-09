#! /usr/bin/env python
# Written on 4/2017 by Sheng-Jun Lin
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import warnings


def polseg_convert(I_map, polI_map, polPA_map,
                   scale_10percent, sampling_interval, seg_color, output_reg):
    """
    PolSeg_convert(I_map, polI_map, polPA_map,
        scale_10percent, sampling_rate, seg_color, output_reg)
    Generate a ds9 region file which contains polarization segments.
    Assume the fits files have 4 axes. (e.g. CASA simulation outputs)
    I_map             [str]: The fits filename of Stokes I.
    polI_map          [str]: The fits filename of Polarized intensity.
    polPA_map         [str]: The fits filename of PA[deg] of polarization segments.
    scale_10percent [float]: Length[arcsec] of 10% polarization segments.
    sampling_interval [int]: (1/sampling rate of segments)[pixel].
    seg_color         [str]: Color of segments.
    output_reg        [str]: Output region filename.
    """
    # Open each fits flie
    I_hdulist = pyfits.open(I_map)
    polI_hdulist = pyfits.open(polI_map)
    polPA_hdulist = pyfits.open(polPA_map)

    # Get the cube and header from HDU lists
    I_data = I_hdulist[0].data
    polI_data = polI_hdulist[0].data
    polI_hd = polI_hdulist[0].header
    polPA_data = polPA_hdulist[0].data

    # CASA simulation outputs have 4 dimension: (Stokes, freq, y, x)
    NS, Nf, Ny, Nx = polI_data.shape
    if NS > 1 or Nf > 1:
        # If the cube is continuum data, NS = 1 (only Stokes I) and Nf = 1
        warnings.warn('n(Stokes)={0}, n(freq)={1}.'.format(
            NS, Nf), RuntimeWarning)
    nx, ny = np.meshgrid(np.arange(Nx), np.arange(Ny))

    # Open the output regoin file and write the header
    reg_file = open(output_reg, "w")
    reg_file.write(
        '# Region file format: DS9 version 4.1\n'
        'global color={0} dashlist=8 3 width=2 font="helvetica 10 normal" '
        'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 '
        'source=1\nfk5\n'.format(seg_color))

    # Get the coordinates (of world coor. system) of polarized intensity map
    wcs = WCS(polI_hd)
    for j in xrange(Ny):
        for i in xrange(Nx):
            # Get the RA and Dec in deg for each pixel
            ra = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)[0]
            dec = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)[1]
            I = I_data[0, 0, j, i]
            polI = polI_data[0, 0, j, i]
            polPA = polPA_data[0, 0, j, i] / 180. * np.pi
            if I > 0. and polI > 0. and not np.isnan(polPA) \
                    and i % sampling_interval == 0 and j % sampling_interval == 0:
                # Calculate the length of segments:
                # 0.5*scale_10percent*10/3600*polI/I*sin(PA) in deg
                dx_half = 0.5 / 360 * scale_10percent * \
                    polI / I * np.sin(polPA) / np.cos(dec / 180. * np.pi)
                dy_half = 0.5 / 360 * scale_10percent * \
                    polI / I * np.cos(polPA)
                reg_file.write("line({0},{1},{2},{3}) # line=0 0\n".format(
                    ra + dx_half, dec + dy_half, ra - dx_half, dec - dy_half))

    reg_file.close()
