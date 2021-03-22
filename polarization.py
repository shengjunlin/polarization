#! /usr/bin/env python
# Modified on 4/2017 by Sheng-Jun Lin from Jia-Wei Wang's script
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import sys, warnings

if (sys.version_info.major >= 3):
    # Python3
    xrange = range

def polseg_convert(I_map, polI_map, polPA_map='',
                   scale_10percent=10., sampling_interval=3,
                   uniform_PA=np.nan, uniform_scale=np.nan, PA_offset=0.,
                   i_hdu=0, i_chan=0, I_clip=0., polI_clip=0.,
                   seg_color='red', output_reg='output.reg', hist_plot=True):
    """
    polseg_convert(I_map, polI_map, polPA_map='',
                   scale_10percent=10., sampling_interval=3,
                   uniform_PA=np.nan, uniform_scale=np.nan, PA_offset=0.,
                   i_hdu=0, i_chan=0, I_clip=0., polI_clip=0.,
                   seg_color='red', output_reg='output.reg', hist_plot=True)

    Use the small-angle approximation to calculate the coordinates of the endpoints of segments,
    and generate a ds9 region file storing the polarization segments.
    The input fits files can have either 4 (e.g., CASA simulation outputs),
                                      or 3 (freq/vel, Dec, RA),
                                      or 2 (Dec, RA) axes.

    Args:
      I_map             [str]: The fits filename of Stokes I map.
      polI_map          [str]: The fits filename of the polarized intensity map.
      polPA_map         [str]: The fits filename of PA [deg] of polarization segments. (Default value = '')
      scale_10percent [float]: The segment length [arcsec] of polarization fraction = 10%. (Default value = 10)
      sampling_interval [int]: (1/sampling rate of segments)[pixel]. (Default value = 3)
      uniform_PA      [float]: Assgin a single PA [deg] to create the segments,
                               which will overwrite "polPA_map" if it isn't a np.nan. (Default value = np.nan)
      uniform_scale   [float]: Create the segments with uniform lengths [arcsec],
                               which will overwrite "scale_10percent" if it isn't a np.nan. (Default value = np.nan)
      PA_offset       [float]: Additional PA offset [deg]. (Default value = 0.)
                               That is, the resulting PA = PA read from polPA_map + PA_offset,
                                                   or PA = uniform_PA + PA_offset.
      i_hdu             [int]: The number in the hdu list. (Default value = 0)
      i_chan            [int]: Channel number. (Default value = 0)
      I_clip          [float]: Exclude pixels in I_map with abs(I) <= I_clip. (Default value = 0)
      polI_clip       [float]: Exclude pixels in polI_map with polI <= polI_clip. (Default value = 0)
      seg_color         [str]: Color of segments. (Defualt value = 'red')
      output_reg        [str]: Output region filename. (Default value = 'output.reg')
      hist_plot        [bool]: Plot a histgram of polarization fraction. (Default value = True)
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

    polper_ls = []
    if np.isnan(uniform_scale):
        # 100% length = (10*scale_10percent[arcsec])/3600 in degree
        segscale_deg = scale_10percent / 360.
    else:
        segscale_deg = uniform_scale / 3600.

    # Get the coordinates (of world coor. system) of polarized intensity map
    wcs = WCS(polI_hd)
    for j in xrange(Ny):
        for i in xrange(Nx):
            # Get the I, polI and polPA for each pixel
            I = I_data[j, i]
            polI = polI_data[j, i]
            if np.isnan(uniform_PA):
                polPA_rad = (polPA_data[j, i] + PA_offset) / 180. * np.pi
            else:
                polPA_rad = (uniform_PA + PA_offset) / 180. * np.pi
            if abs(I) > I_clip and polI > polI_clip and not np.isnan(polPA_rad) \
                    and i % sampling_interval == 0 and j % sampling_interval == 0:
                # Get the RA and Dec in deg for each pixel
                try:
                    # 4-dim
                    ra_deg, dec_deg, _, _ = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0, 0)
                except:
                    try:  # 3-dim
                        ra_deg, dec_deg, _ = wcs.all_pix2world(nx[j, i], ny[j, i], 0, 0)
                    except:  # 2-dim
                        ra_deg, dec_deg = wcs.all_pix2world(nx[j, i], ny[j, i], 0)
                # Calculate an half of length of segments:
                polper = polI/I
                polper_ls.append(polper)
                # 0.5 * (polI/I) * sin[or cos](PA) * (10*scale_10percent[arcs])/3600 is in deg
                # dRA is corrected by cos(Dec)
                dx_half_deg = 0.5 * segscale_deg * polper * np.sin(polPA_rad) / np.cos(dec_deg / 180. * np.pi)
                dy_half_deg = 0.5 * segscale_deg * polper * np.cos(polPA_rad)
                reg_file.write("line({0},{1},{2},{3}) # line=0 0\n".format(
                    ra_deg + dx_half_deg, dec_deg + dy_half_deg, ra_deg - dx_half_deg, dec_deg - dy_half_deg))

    reg_file.close()

    if hist_plot:
        from matplotlib.pyplot import hist, xlabel, ylabel, show
        hist(polper_ls)
        xlabel('polarized percentage')
        ylabel('N')
        show()
