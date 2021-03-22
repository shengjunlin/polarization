#! /usr/bin/env python
# Modified on 4/2017 by Sheng-Jun Lin from Jia-Wei Wang's script
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import warnings


def polseg_convert(I_map, polI_map, polPA_map, scale_10percent, i_hdu=0, i_chan=0, sampling_interval=3,
                   I_clip=0., polI_clip=0., seg_color='red', uniformPA=np.nan, output_reg='output.reg',
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
    i_hdu             [int]: The number in the hdu list. Default: 0.
    i_chan            [int]: Channel number. Default: 0.
    sampling_interval [int]: (1/sampling rate of segments)[pixel]. Default: 3.
    I_clip          [float]: Exclude pixels in I_map with values <= I_clip. Default: 0.
    polI_clip       [float]: Exclude pixels in polI_map with values <= polI_clip. Default: 0.
    seg_color         [str]: Color of segments. Defualt: 'red'.
    uniformPA
    output_reg        [str]: Output region filename. Default: 'output.reg'.
    hist_plot     [boolean]: Plot a histgram of pol. percentage. Default: True.
    """
    # Open each fits flie
    I_hdulist = pyfits.open(I_map)
    polI_hdulist = pyfits.open(polI_map)
    polPA_hdulist = pyfits.open(polPA_map)

    # Get the cube and header from HDU lists
    I_data = I_hdulist[i_hdu].data
    polI_data = polI_hdulist[i_hdu].data
    polI_hd = polI_hdulist[i_hdu].header
    polPA_data = polPA_hdulist[i_hdu].data

    try:
        # CASA simulation outputs have 4 dimension: (Stokes, freq, y, x)
        NS, Nf, Ny, Nx = polI_data.shape
        if NS > 1 or Nf > 1:
            # If the cube is continuum data, NS = 1 (only Stokes I) and Nf = 1
            warnings.warn('n(Stokes)={0}, n(freq)={1}.'.format(
                NS, Nf), RuntimeWarning)
        I_data = np.squeeze(I_data[0, i_chan, :, :])
        polI_data = np.squeeze(polI_data[0, i_chan, :, :])
        polPA_data = np.squeeze(polPA_data[0, i_chan, :, :])
    except:
        try:
            # 3-dim fits
            Nf, Ny, Nx = polI_data.shape
            I_data = np.squeeze(I_data[i_chan, :, :])
            polI_data = np.squeeze(polI_data[i_chan, :, :])
            polPA_data = np.squeeze(polPA_data[i_chan, :, :])
        except:
            # 2-dim fits
            Ny, Nx = polI_data.shape

    nx, ny = np.meshgrid(np.arange(Nx), np.arange(Ny))

    # Open the output regoin file and write the header
    reg_file = open(output_reg, "w")
    reg_file.write(
        '# Region file format: DS9 version 4.1\n'
        'global color={0} dashlist=8 3 width=1 font="helvetica 10 normal" '
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
            I = I_data[j, i]
            polI = polI_data[j, i]
            if np.isnan(uniformPA):
                polPA_rad = polPA_data[j, i] / 180. * np.pi
            else:
                polPA_rad = uniformPA / 180. * np.pi
            if abs(I) > I_clip and polI > polI_clip and not np.isnan(polPA_rad) \
                    and i % sampling_interval == 0 and j % sampling_interval == 0:
                # Get the RA and Dec in deg for each pixel
                try:
                    # 4-dim
                    ra_deg, dec_deg, _, _ = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)
                except:
                    try:
                        ra_deg, dec_deg, _ = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0)
                    except:
                        ra_deg, dec_deg = wcs.all_pix2world(nx[j, i], ny[j, i], 0)
                # ra_deg = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0)[0]
                # dec_deg = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0)[1]
                # Calculate an half of length of segments:
                polper = polI/I
                polper_ls.append(polper)
                # 0.5 * (polI/I) * sin[or cos](PA) * (10*scale_10percent[arcs])/3600 is in deg
                # dRA is corrected by cos(Dec)
                dx_half_deg = 0.5 * polper * np.sin(polPA_rad) / np.cos(dec_deg / 180. * np.pi) \
                    * scale_10percent / 360.
                dy_half_deg = 0.5 * polper * np.cos(polPA_rad) \
                    * scale_10percent / 360.
                reg_file.write("line({0},{1},{2},{3}) # line=0 0\n".format(
                    ra_deg + dx_half_deg, dec_deg + dy_half_deg, ra_deg - dx_half_deg, dec_deg - dy_half_deg))

    reg_file.close()

    if hist_plot:
        hist(polper_ls)
        xlabel('polarized percentage')
        ylabel('N')
        show()
