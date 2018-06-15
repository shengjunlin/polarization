#! /usr/bin/env python
# Modified on 4/2017 by Sheng-Jun Lin from Jia-Wei Wang's script
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import warnings


def polseg_convert(I_map, polI_map, polPA_map, scale_10percent, sampling_interval=3,
                   I_clip=0., polI_clip=0., seg_color='r', output_reg='output.reg',
                   hist_plot=True):
    """
    polseg_convert(I_map, polI_map, polPA_map,
        scale_10percent, sampling_interval, seg_color, output_reg)
    Generate a ds9 region file which contains polarization segments.
    Assume the fits files have 4 axes. (e.g. CASA simulation outputs)
    I_map             [str]: The fits filename of Stokes I.
    polI_map          [str]: The fits filename of Polarized intensity.
    polPA_map         [str]: The fits filename of PA[deg] of polarization segments.
    scale_10percent [float]: Length[arcsec] of 10% polarization segments.
    sampling_interval [int]: (1/sampling rate of segments)[pixel]. Default: 3.
    I_clip          [float]: Exclude pixels in I_map with values <= I_clip. Default: 0.
    polI_clip       [float]: Exclude pixels in polI_map with values <= polI_clip. Default: 0.
    seg_color         [str]: Color of segments. Defualt: 'red'.
    output_reg        [str]: Output region filename. Default: 'output.reg'.
    hist_plot     [boolean]: Plot a histgram of pol. percentage. Default: True.
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

    if hist_plot:
        from matplotlib.pyplot import hist, xlabel, ylabel, show
    polper_ls = []

    # Get the coordinates (of world coor. system) of polarized intensity map
    wcs = WCS(polI_hd)
    for j in xrange(Ny):
        for i in xrange(Nx):
            # Get the I, polI and polPA for each pixel
            I = I_data[0, 0, j, i]
            polI = polI_data[0, 0, j, i]
            polPA = polPA_data[0, 0, j, i] / 180. * np.pi
            if I > I_clip and polI > polI_clip and not np.isnan(polPA) \
                    and i % sampling_interval == 0 and j % sampling_interval == 0:
                # Get the RA and Dec in deg for each pixel
                ra = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)[0]
                dec = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)[1]
                # Calculate an half of length of segments:
                polper = polI/I
                polper_ls.append(polper)
                # 0.5 * (polI/I) * sin[or cos](PA) * (10*scale_10percent[arcs])/3600 in deg
                dx_half = 0.5 / 360 * scale_10percent * \
                    polper * np.sin(polPA) / np.cos(dec / 180. * np.pi) # dRA is corrected by cos(Dec)
                dy_half = 0.5 / 360 * scale_10percent * \
                    polper * np.cos(polPA)
                reg_file.write("line({0},{1},{2},{3}) # line=0 0\n".format(
                    ra + dx_half, dec + dy_half, ra - dx_half, dec - dy_half))

    reg_file.close()

    if hist_plot:
        hist(polper_ls)
        xlabel('polarized percentage')
        ylabel('N')
        show()
