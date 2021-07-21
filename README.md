# polarization
A Python function related to polarization.

Description : First modify Jia-Wei Wang's script that reading
I, Q, U images to generate a ds9 region file of polarization
segments. Then extend it to a few functions reading not just
ALMA CASA images and also JCMT SCUBA2-POL2 catalogues.

Date :

2017-04 Jia-Wei's original script.    
2021-07 POL2 catalogue extension.    

In the python shell, import polarization with the following code:  

    import sys
    sys.path.append('the path of the dir containing the script')
    from polarization import *


def write_regfilehd(reg_file, seg_color='black', seg_width=1):  

    """Write the standard header of the ds9 region file object.

    Parameters
    ----------
    reg_file : file object
        Ds9 region file object.
    seg_color : str
        Segment color. (Default value = 'black')
    seg_width : int
        Segment width. (Default value = 1)

    Returns
    -------
    reg_file : file object
        Ds9 region file object.

    """


def polseg_to_regfile(reg_file, RA, Dec, RADec_unit='deg', p=None, ten_percent_scale_asec=None, uniform_scale_asec=None, PA_deg=None, uniform_PA_deg=None, PA_offset=0.):  

    """Write polarisation segments into the ds9 region file object.
    This is modified by Jia-Wei Wang's script.

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


def polseg_from_POL2_cat(reg_filename, POL2_cat_filename='', data=None, mask=None, ten_percent_scale_asec=None, uniform_scale_asec=None, uniform_PA_deg=None, PA_offset=90., seg_color='black', seg_width=1):  

    """Convert a POL2 catlogue into a segment ds9 region file with the given mask.
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


    """

def polseg_from_2darray(reg_filename, header, I_data, PI_data, PA_data_deg=None, sampling_interval_px=3, ten_percent_scale_asec=None, uniform_scale_asec=None, uniform_PA_deg=None, PA_offset=0., I_clip=0., PI_clip=0., mask_func=None, seg_color='black', seg_width=1):  

    """Convert the given I/PI/PA arrays into a segment ds9 region file with the given mask_func.
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
        set by "ten_percent_scale_asec". (Default value = None)
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

    Returns
    -------
    reg_file : file object
        Ds9 region file object.
    ls_p : list
        The list of the polarisation fractions.

    """


def polseg_from_images(reg_filename, I_map, PI_map, PA_map_deg='', i_hdu=0, i_chan=0, sampling_interval_px=3, ten_percent_scale_asec=None, uniform_scale_asec=None, uniform_PA_deg=None, PA_offset=0., I_clip=0., PI_clip=0., mask_func=None, seg_color='black', seg_width=1, hist_plot=True):  

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
        set by "ten_percent_scale_asec". (Default value = None)
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
    hist_plot : bool
        Plot a histgram of polarization fraction. (Default value = True)

    """


def polseg_convert(I_map, polI_map, polPA_map=None, scale_10percent=10., sampling_interval=3, uniform_PA=None, uniform_scale=None, PA_offset=0., i_hdu=0, i_chan=0, I_clip=0., polI_clip=0., seg_color='red', output_reg='output.reg', hist_plot=True):  

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

