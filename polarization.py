#! /usr/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin
# Email : shengjunlin@asiaa.sinica.edu.tw
# Description : First modify Jia-Wei Wang's script that reading
# I, Q, U images to generate a ds9 region file of polarization
# segments. Then extend it to a few functions reading not just
# ALMA CASA images and also JCMT SCUBA2-POL2 catalogues.
# Date : 2017-04 Jia-Wei's original script.
#        2021-07 POL2 catalogue extension. Restructure functions.
# #########################################################
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import warnings
import sys

if (sys.version_info.major >= 3):
    # Python3
    xrange = range

def get_dir_filename(filepath):

    """Return the directories and filename.
    Note: If input a path of mutiple directories (end with '/'),
          get_dir_filename returns (rest_dirs_name, last_dirname).

    Args:
      filepath: string.
    Returns:
      dirs_name: string.
      filename: string.
    """
    i = filepath.rfind('/', 0, -1) + 1  # Choose -1 to avoid '/' at the end
    dirs_name = filepath[:i]
    filename = filepath[i:]
    return (dirs_name, filename)


def write_regfilehd(reg_file, seg_color='black', seg_width=1, frame='fk5'):
    """Write the standard header of the ds9 region file object.

    Parameters
    ----------
    reg_file : file object
        Ds9 region file object.
    seg_color : str
        Segment color. (Default value = 'black')
    seg_width : int
        Segment width. (Default value = 1)
    frame : {'fk5', 'icrs', 'fk4'}
        Please choose 'icrs' but not 'fk5' for ALMA data.
        (Default value = 'fk5')

    Returns
    -------
    reg_file : file object
        Ds9 region file object.

    """

    reg_file.write(
        '# Region file format: DS9 version 4.1\n'
        'global color={0} dashlist=8 3 width={1} font="helvetica 10 normal" '
        'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 '
        'source=1\n{2}\n'.format(seg_color, seg_width, frame))
    return reg_file


def polseg_to_regfile(reg_file, RA, Dec, RADec_unit='deg',
        p=None, ten_percent_scale_asec=None, uniform_scale_asec=None,
        PA_deg=None, uniform_PA_deg=None, PA_offset=0.):
    """Write polarisation segments into the ds9 region file object.
    This is modified from Jia-Wei Wang's script.

    Parameters
    ----------
    reg_file : file object
        Ds9 region file object.
    RA, Dec : float
        The coordinates in either deg or rad.
    RADec_unit : {'deg', 'rad'}
        The units of RA and Dec. (Default value = 'deg')
    p : float
        Polarisation fraction, PI/I. (Default value = None)
    ten_percent_scale_asec : float
        The segment length in arcsec at p = 10%. (Default value = None)
    uniform_scale_asec : float
        Assign uniform-length segments, overwriting the p-scaled segments
        set by "p" and "ten_percent_scale_asec". (Default value = None)
    PA_deg : float
        Polarisation angle in degree. (Default value = None)
    uniform_PA_deg : float
        Assign an uniform PA for plotting, overwriting "PA_deg".
        (Default value = None)
    PA_offset : float
        Additional offset in deg to "PA_deg" or "uniform_PA_deg".
        (Default value = 0.)

    Returns
    -------
    reg_file : file object
        Ds9 region file object.

    """

    # Determine the segment length in deg
    if uniform_scale_asec is None:
        seg_length_deg = p * (10 * ten_percent_scale_asec / 3600.)
    else:
        seg_length_deg = uniform_scale_asec / 3600.

    # Determine PA in deg
    if uniform_PA_deg is None:
        PA_rad = (PA_deg + PA_offset) / 180. * np.pi
    else:
        PA_rad = (uniform_PA_deg + PA_offset) / 180. * np.pi

    # RA and Dec in rad or deg
    if RADec_unit == 'deg':
        RA_deg = RA
        Dec_deg = Dec
        Dec_rad = Dec / 180. * np.pi
    elif RADec_unit == 'rad':
        RA_deg = RA / np.pi * 180.
        Dec_deg = Dec / np.pi * 180.
        Dec_rad = Dec
    else:
        raise RuntimeError('polseg_to_regfile: The units of RA and Dec are unkown.')

    # Calculate an half of length of segments:
    # 0.5 * seg_length_deg * [sin|cos](PA), and dRA is corrected by cos(Dec)
    half_dx_deg = 0.5 * seg_length_deg * np.sin(PA_rad) / np.cos(Dec_rad)
    half_dy_deg = 0.5 * seg_length_deg * np.cos(PA_rad)
    reg_file.write("line({0},{1},{2},{3}) # line=0 0\n".format(
        RA_deg  + half_dx_deg,
        Dec_deg + half_dy_deg,
        RA_deg  - half_dx_deg,
        Dec_deg - half_dy_deg))
    return reg_file


def polseg_from_POL2_cat(reg_filename, POL2_cat_filename='', data=None, mask=None,
        ten_percent_scale_asec=None, uniform_scale_asec=None,
        uniform_PA_deg=None, PA_offset=90., seg_color='black', seg_width=1, frame='fk5'):
    """Convert a POL2 catlogue into a ds9 region file with the given mask.
    The endpoints of segments are calculated with the small-angle approximation.
    To obtain the B-field segments, the default of "PA_offset" sets to 90 degree.

    Parameters
    ----------
    reg_filename : str
        Output ds9 region filename.
    POL2_cat_filename : str
        POL2 catalgue filename. (Default value = '')
    data : BinTableHDU object, optional
        Assign a BinTableHDU object instead of "POL2_cat_filename" by
            data = fits.open(POL2cat_filename, memmap=True)[1].data
        (Default value = None)
    mask : bool 1d-array
        The mask to filter the POL2 catalogue. The lenths of "data" and
        "mask" are equal. There is only one bulit-in filter: I >= 0.
        (Default value = None)
    ten_percent_scale_asec : float
        The segment length in arcsec at p = 10%. (Default value = None)
    uniform_scale_asec : float
        Assign uniform-length segments, overwriting the p-scaled segments
        set by "ten_percent_scale_asec". (Default value = None)
    uniform_PA_deg : float
        Assign an uniform PA for plotting, overwriting PA in the catlogue.
        (Default value = None)
    PA_offset : float
        Additional offset in deg to PA in the catlogue or "uniform_PA_deg".
        (Default value = 90.)
    seg_color : str
        Segment color. (Default value = 'black')
    seg_width : int
        Segment width. (Default value = 1)
    frame : {'fk5', 'icrs', 'fk4'}
        Please choose 'icrs' but not 'fk5' for ALMA data.
        (Default value = 'fk5')


    """
    # from astropy.table import Table
    # Note: Not sure what is the advantage to use Table?

    if data is None:
        data = fits.open(POL2_cat_filename, memmap=True)[1].data

    if mask is None:
        table = data
        # table = Table(data)
    else:
        table = data[mask]
        # table = Table(data[mask])

    table['P'] /= 100.  # POL2 catalogue stores 100*P
    table['DP'] /= 100.  # POL2 catalogue stores 100*DP

    # Open the output regoin file and write the header
    with open(reg_filename, "w") as reg_file:
        reg_file = write_regfilehd(reg_file, seg_color, seg_width, frame)

        for i in xrange(len(table)):
            RA_rad = table['RA'][i]  # rad
            Dec_rad = table['DEC'][i]  # rad
            p = table['P'][i]
            dp = table['DP'][i]
            PA_deg = table['ANG'][i]  # deg
            I = table['I'][i]  # mJy/bm

            if I >= 0.:
                reg_file = \
                    polseg_to_regfile(reg_file,
                                      RA_rad,
                                      Dec_rad,
                                      RADec_unit='rad',
                                      p=p,
                                      ten_percent_scale_asec=ten_percent_scale_asec,
                                      uniform_scale_asec=uniform_scale_asec,
                                      PA_deg=PA_deg,
                                      uniform_PA_deg=uniform_PA_deg,
                                      PA_offset=PA_offset)


def polseg_from_2darray(reg_filename, header,
        I_data=None, PI_data=None, PA_data_deg=None, sampling_interval_px=3,
        ten_percent_scale_asec=None, uniform_scale_asec=None,
        uniform_PA_deg=None, PA_offset=0.,
        I_clip=0., PI_clip=0., mask_func=None, seg_color='black', seg_width=1, frame='fk5'):
    """Convert the given I/PI/PA arrays into a ds9 region file with the given mask_func.
    The endpoints of segments are calculated with the small-angle approximation.

    Parameters
    ----------
    reg_filename : str
        Output ds9 region filename.
    header : header of an HDU
        A header for reading the wcs. Assume that all of the arrays share the same one.
    I_data : 2d-array
        The Stokes I map.
    PI_data : 2d-array
        The polarised intensity map. I_data and PI_data should have the same units.
    PA_data_deg : 2d-array
        The polarisation angle map in deg. (Default value = None)
    sampling_interval_px : int
        The sampling interval in pixel for segments ( = 1/sampling rate ).
        (Default value = 3)
    ten_percent_scale_asec : float
        The segment length in arcsec at p = 10%. (Default value = None)
    uniform_scale_asec : float
        Assign uniform-length segments, overwriting the p-scaled segments
        set by "ten_percent_scale_asec". "I_data" and "PI_data" are useless,
        and "I_clip" and "PI_clip" will set to -np.inf.
        (Default value = None)
    uniform_PA_deg : float
        Assign an uniform PA for plotting, overwriting PA from "PA_data_deg".
        (Default value = None)
    PA_offset : float
        Additional offset in deg to "PA_data_deg" or "uniform_PA_deg".
        (Default value = 0.)
    I_clip : float
        Exclude pixels in I_data with abs(I) < I_clip. (Default value = 0)
    PI_clip : float
        Exclude pixels in PI_data with PI < PI_clip. (Default value = 0)
    mask_func : callable
        A True/False-valued function of multiple scalar variables:
        i, j, RA_deg, Dec_deg, I, PI, PA_deg, where i and j are python-like
        indices of x (RA) and y (Dec) positions.
        So explicitly define it like:
            def mask_func(i, j, RA_deg, Dec_deg, I, PI, PA_deg):
                if ... :
                    return True
                else:
                    return False
    seg_color : str
        Color of segments. (Defualt value = 'black')
    seg_width : int
        Segment width. (Default value = 1)
    frame : {'fk5', 'icrs', 'fk4'}
        Please choose 'icrs' but not 'fk5' for ALMA data.
        (Default value = 'fk5')

    Returns
    -------
    reg_file : file object
        Ds9 region file object.
    ls_p : list
        The list of the polarisation fractions.

    """
    # Get the coordinates from the input header
    wcs = WCS(header)
    ndim = header['NAXIS']
    if uniform_scale_asec is None:
        Ny, Nx = I_data.shape
    else:
        Ny, Nx = PA_data_deg.shape
        I_clip = -np.inf
        PI_clip = -np.inf
    nx, ny = np.meshgrid(np.arange(Nx), np.arange(Ny))

    # Get the RA and Dec in deg for each pixel
    if ndim == 4:
        coor_list = wcs.all_pix2world(nx, ny, 0, 0, 0)
        # coor_list is a 4-element list.
        # coor_list = [RA_deg_2Dmap, Dec_deg_2Dmap, Freq_Hz_2Dmap, Stokes_2Dmap]
    elif ndim == 3:
        coor_list = wcs.all_pix2world(nx, ny, 0, 0)
        # coor_list is a 3-element list.
        # coor_list = [RA_deg_2Dmap, Dec_deg_2Dmap, Freq_Hz_2Dmap OR Velo_kms_2Dmap]
    elif ndim == 2:
        coor_list = wcs.all_pix2world(nx, ny, 0)
        # coor_list is a 2-element list.
        # coor_list = [RA_deg_2Dmap, Dec_deg_2Dmap]
    else:
        raise RuntimeError('polseg_from_2darray: '
                'The dimension of the input header is weird.')

    # Each element in "coor_list" is a 2d-array with the dim of Nx by Ny.
    RA_data_deg  = coor_list[0]  # deg
    Dec_data_deg = coor_list[1]  # deg

    # If the mask function is undefined
    if mask_func is None:
        mask_func = lambda i, j, RA_deg, Dec_deg, I, PI, PA_deg: True

    # Open the output regoin file and write the header
    with open(reg_filename, "w") as reg_file:
        reg_file = write_regfilehd(reg_file, seg_color, seg_width, frame)
        ls_p = []
        for j in xrange(Ny):
            for i in xrange(Nx):
                # Get the values for each pixel
                RA_deg = RA_data_deg[j, i]
                Dec_deg = Dec_data_deg[j, i]
                if uniform_scale_asec is None:
                    I = I_data[j, i]
                    PI = PI_data[j, i]
                    p = PI / I
                else:
                    I = np.inf
                    PI = np.inf
                    p = None
                PA_deg = PA_data_deg[j, i]
                if abs(I) >= I_clip and \
                        PI >= PI_clip and \
                    (not np.isnan(PA_deg)) and \
                    mask_func(i, j, RA_deg, Dec_deg, I, PI, PA_deg) and \
                    i % sampling_interval_px == 0 and \
                    j % sampling_interval_px == 0:

                    reg_file = \
                        polseg_to_regfile(reg_file,
                                          RA_deg,
                                          Dec_deg,
                                          RADec_unit='deg',
                                          p=p,
                                          ten_percent_scale_asec=ten_percent_scale_asec,
                                          uniform_scale_asec=uniform_scale_asec,
                                          PA_deg=PA_deg,
                                          uniform_PA_deg=uniform_PA_deg,
                                          PA_offset=PA_offset)
                    ls_p.append(p)

    return (reg_file, ls_p)


def polseg_from_images(reg_filename,
        I_map=None, PI_map=None, PA_map_deg='',
        i_hdu=0, i_chan=0,
        sampling_interval_px=3,
        ten_percent_scale_asec=None, uniform_scale_asec=None,
        uniform_PA_deg=None, PA_offset=0.,
        I_clip=0., PI_clip=0., mask_func=None,
        seg_color='black', seg_width=1, frame='icrs',
        hist_plot=True):

    """Use the small-angle approximation to calculate the coordinates of
    the endpoints of segments, and generate a ds9 region file storing
    the polarization segments.
    The input fits files can have either 4 (e.g., CASA simulation outputs),
                                      or 3 (freq/vel, Dec, RA),
                                      or 2 (Dec, RA) axes.

    Parameters
    ----------
    reg_filename : str
        Output ds9 region filename.
    I_map : str
        The fits filename of Stokes I map.
    PI_map : str
        The fits filename of the polarised intensity map.
    PA_map_deg : str
        The fits filename of the PA map in deg. (Default value = '')
    i_hdu : int
        The number in the hdu list. (Default value = 0)
    i_chan : int
        Channel number. (Default value = 0)
    sampling_interval_px : int
        The sampling interval in pixel for segments ( = 1/sampling rate ).
        (Default value = 3)
    ten_percent_scale_asec : float
        The segment length in arcsec at p = 10%. (Default value = None)
    uniform_scale_asec : float
        Assign uniform-length segments, overwriting the p-scaled segments
        set by "ten_percent_scale_asec". "I_map" and "PI_map" are useless.
        (Default value = None)
    uniform_PA_deg : float
        Assign an uniform PA for plotting, overwriting PA from "PA_map_deg".
        (Default value = None)
    PA_offset : float
        Additional offset in deg to "PA_map_deg" or "uniform_PA_deg".
        (Default value = 0.)
    I_clip : float
        Exclude pixels in I_map with abs(I) < I_clip. (Default value = 0)
    PI_clip : float
        Exclude pixels in PI_map with PI < PI_clip. (Default value = 0)
    mask_func : callable
        A True/False-valued function of multiple scalar variables:
        i, j, RA_deg, Dec_deg, I, PI, PA_deg, where i and j are python-like
        indices of x (RA) and y (Dec) positions.
        So explicitly define it like:
            def mask_func(i, j, RA_deg, Dec_deg, I, PI, PA_deg):
                if ... :
                    return True
                else:
                    return False
    seg_color : str
        Color of segments. (Defualt value = 'black')
    seg_width : int
        Segment width. (Default value = 1)
    frame : {'fk5', 'icrs', 'fk4'}
        Please choose 'icrs' but not 'fk5' for ALMA data.
        (Default value = 'icrs')
    hist_plot : bool
        Plot a histgram of polarization fraction. (Default value = True)

    """
    # Open each fits flie
    # Get the cube and header from HDU lists
    if uniform_PA_deg is None:
        PA_hdulist = fits.open(PA_map_deg)
        PA_data_deg = PA_hdulist[i_hdu].data
    else:
        PA_data_deg = None

    if uniform_scale_asec is None:
        I_hdulist = fits.open(I_map)
        PI_hdulist = fits.open(PI_map)
        hd = I_hdulist[i_hdu].header
        I_data = I_hdulist[i_hdu].data
        PI_data = PI_hdulist[i_hdu].data
        ndim = len(I_data.shape)
    else:  # PA map must be given
        hd = PA_hdulist[i_hdu].header
        ndim = len(PA_data_deg.shape)

    if ndim == 4:
        # CASA simulation outputs have 4 dimension: (Stokes, freq, y, x)
        if uniform_scale_asec is None:
            NS, Nf, Ny, Nx = I_data.shape
        else:
            NS, Nf, Ny, Nx = PA_data_deg.shape
        if NS > 1 or Nf > 1:
            # If the cube is continuum data, NS = 1 (only Stokes I) and Nf = 1
            warnings.warn('polseg_from_images: n(Stokes)={0}, n(freq)={1}.'.format(
                NS, Nf), RuntimeWarning)
        if uniform_scale_asec is None:
            I_data = np.squeeze(I_data[0, i_chan, :, :])
            PI_data = np.squeeze(PI_data[0, i_chan, :, :])
        else:
            I_data = None
            PI_data = None
        if uniform_PA_deg is None:
            PA_data_deg = np.squeeze(PA_data_deg[0, i_chan, :, :])
    elif ndim == 3:
        # 3-dim fits
        if uniform_scale_asec is None:
            I_data = np.squeeze(I_data[i_chan, :, :])
            PI_data = np.squeeze(PI_data[i_chan, :, :])
        else:
            I_data = None
            PI_data = None
        if uniform_PA_deg is None:
            PA_data_deg = np.squeeze(PA_data_deg[i_chan, :, :])
    elif ndim == 2:
        # 2-dim fits
        pass
    else:
        raise RuntimeError('polseg_from_images: '
                'The dimension of one of the input fits files is worng.')

    reg_file, ls_p = \
        polseg_from_2darray(reg_filename,
                            hd,
                            I_data,
                            PI_data,
                            PA_data_deg,
                            sampling_interval_px,
                            ten_percent_scale_asec,
                            uniform_scale_asec,
                            uniform_PA_deg,
                            PA_offset,
                            I_clip,
                            PI_clip,
                            mask_func,
                            seg_color,
                            seg_width,
                            frame)

    if uniform_scale_asec is None:
        I_hdulist.close()
        PI_hdulist.close()
    if uniform_PA_deg is None:
        PA_hdulist.close()

    if hist_plot:
        from matplotlib.pyplot import hist, xlabel, ylabel, show
        hist(ls_p)
        xlabel('polarized fraction')
        ylabel('N')
        show()


def pol2mapCat_to_fitsImages(POL2_cat_filepath,
        template_fits_filepath,
        output_dirpath='./',
        jy=True,
        index_offset=0):
    """Load the cat generated by pol2map, and use the header
    from the template fits file to create fits files for
    I, DI, Q, DQ, U, DU, P, DP, ANG, DANG, PI, DPI maps,
    where P, DP, PI, DPI maps are denpending on the debias mode.

    The output fits filename is
        `template_fits_filename`_[I|DI|Q|DQ|U|DU|ANG|DANG].fits,
    and
        `POL2_cat_filename`_[P|DP|PI|DPI].fits

    Parameters
    ----------
    POL2_cat_filepath : str

    template_fits_filepath : str

    output_dirpath : str
        (Default value = './')
    jy : bool
        The units of I, DI, Q, DQ, U, DU, PI, DPI are mJy/beam if True.
        Otherwise, their units are pW in `POL2_cat_filepath`.
        (Default value = True)
    index_offset : int
        A testing parameter. L1512 seems 0, L1498 seems 1
    """

    if not jy:
        print('Error!')

    with fits.open(POL2_cat_filepath, memmap=True) as cat_hdul:
        data = cat_hdul[1].data
        X = data['X']
        length = len(X)
        Y = data['Y']
        # X, Y = 0.5, 0.5 is crpix_x, crpix_y
        I = data['I'] / 1e3  # Jy/bm
        DI = data['DI'] / 1e3  # Jy/bm
        Q = data['Q'] / 1e3  # Jy/bm
        DQ = data['DQ'] / 1e3  # Jy/bm
        U = data['U'] / 1e3  # Jy/bm
        DU = data['DU'] / 1e3  # Jy/bm
        PI = data['PI'] / 1e3  # Jy/bm
        DPI = data['DPI'] / 1e3  # Jy/bm
        p = data['P'] / 100.
        dp = data['DP'] / 100.
        ANG = data['ANG']
        DANG = data['DANG']
        pSNR = p/dp

    with fits.open(template_fits_filepath) as fits_hdul:
        hd = fits_hdul[0].header
        crpix_x = hd['CRPIX1']
        crpix_y = hd['CRPIX2']

        data = fits_hdul[0].data

        # I
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = I[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_I.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # DI
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = DI[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_DI.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # Q
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = Q[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_Q.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # DQ
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = DQ[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_DQ.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # U
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = U[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_U.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # DU
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = DU[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_DU.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # PI
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = PI[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_PI.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # DPI
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = DPI[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_DPI.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # p
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = p[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_p.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # dp
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = dp[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_dp.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # PA
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = ANG[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_PA.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

        # DPA
        data *= 0.0
        for i in range(length):
            data[0, int(crpix_y+Y[i]-0.5-index_offset), int(crpix_x+X[i]-0.5-index_offset)] = DANG[i]
        _, filename = get_dir_filename(template_fits_filepath.replace('.fits', '_DPA.fits'))
        filepath = output_dirpath + filename
        fits_hdul.writeto(filepath)

def polseg_convert(I_map, polI_map, polPA_map=None,
                   scale_10percent=10., sampling_interval=3,
                   uniform_PA=None, uniform_scale=None, PA_offset=0.,
                   i_hdu=0, i_chan=0, I_clip=0., polI_clip=0.,
                   seg_color='red', output_reg='output.reg', hist_plot=True):

    """A wrapper for "polseg_from_images" to be compatiable with the old scripts
    that use the function, "polseg_convert".

    Parameters
    ----------
    I_map : str
        The fits filename of Stokes I map.
    polI_map : str
        The fits filename of the polarised intensity map.
    polPA_map : str
        The fits filename of the PA map in deg. (Default value = '')
    scale_10percent : float
        The segment length in arcsec at p = 10%. (Default value = None)
    sampling_interval : int
        The sampling interval in pixel for segments ( = 1/sampling rate ).
        (Default value = 3)
    uniform_PA : float
        Assign an uniform PA for plotting, overwriting PA from "polPA_map".
        (Default value = None)
    uniform_scale : float
        Assign uniform-length segments, overwriting the p-scaled segments
        set by "scale_10percent". (Default value = None)
    PA_offset : float
        Additional offset in deg to "PA_map_deg" or "uniform_PA_deg".
        (Default value = 0.)
    i_hdu : int
        The number in the hdu list. (Default value = 0)
    i_chan : int
        Channel number. (Default value = 0)
    I_clip : float
        Exclude pixels in I_map with abs(I) < I_clip. (Default value = 0)
    polI_clip : float
        Exclude pixels in polI_map with polI < polI_clip. (Default value = 0)
    seg_color : str
        Color of segments. (Defualt value = 'red')
    output_reg : str
        Output ds9 region filename.
    hist_plot : bool
        Plot a histgram of polarization fraction. (Default value = True)

    """
    polseg_from_images(reg_filename=output_reg,
        I_map=I_map,
        PI_map=polI_map,
        PA_map_deg=polPA_map,
        i_hdu=i_hdu,
        i_chan=i_chan,
        sampling_interval_px=sampling_interval,
        ten_percent_scale_asec=scale_10percent,
        uniform_scale_asec=uniform_scale,
        uniform_PA_deg=uniform_PA,
        PA_offset=PA_offset,
        I_clip=I_clip,
        PI_clip=polI_clip,
        mask_func=None,
        seg_color=seg_color,
        seg_width=1,
        hist_plot=hist_plot)

