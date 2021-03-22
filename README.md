# polarization
A Python function related to polarization.

In the python shell, import polarization with the following code:  

    import sys
    sys.path.append('the path of the dir containing the script')
    from polarization import *

**polseg_convert()**

Syntax:

    polseg_convert(I_map, polI_map, polPA_map='',
                   scale_10percent=10., sampling_interval=3,
                   uniform_PA=np.nan, uniform_scale=np.nan, PA_offset=0.,
                   i_hdu=0, i_chan=0, I_clip=0., polI_clip=0.,
                   seg_color='red', output_reg='output.reg', hist_plot=True)

Arguments:

I_map             [str]: The fits filename of Stokes I map.    
polI_map          [str]: The fits filename of the polarized intensity map.    
polPA_map         [str]: The fits filename of PA [deg] of polarization segments. (Default value = '')    
scale_10percent [float]: The segment length [arcsec] of polarization fraction = 10%. (Default value = 10)    
sampling_interval [int]: (1/sampling rate of segments)[pixel]. (Default value = 3)    
uniform_PA      [float]: Assgin a single PA [deg] to create the segments, which will overwrite "polPA_map" if it isn't a np.nan. (Default value = np.nan)    
uniform_scale   [float]: Create the segments with uniform lengths [arcsec], which will overwrite "scale_10percent" if it isn't a np.nan. (Default value = np.nan)    
PA_offset       [float]: Additional PA offset [deg]. That is, the resulting PA = PA read from polPA_map + PA_offset, or PA = uniform_PA + PA_offset. (Default value = 0.)    
i_hdu             [int]: The number in the hdu list. (Default value = 0)    
i_chan            [int]: Channel number. (Default value = 0)    
I_clip          [float]: Exclude pixels in I_map with abs(I) <= I_clip. (Default value = 0)    
polI_clip       [float]: Exclude pixels in polI_map with polI <= polI_clip. (Default value = 0)    
seg_color         [str]: Color of segments. (Defualt value = 'red')    
output_reg        [str]: Output region filename. (Default value = 'output.reg')    
hist_plot        [bool]: Plot a histgram of polarization fraction. (Default value = True)    

Use the small-angle approximation to calculate the coordinates of the endpoints of segments, and generate a ds9 region file storing the polarization segments.

The input fits files can have either 4 (e.g., CASA simulation outputs), 3 (freq/vel, Dec, RA), or 2 (Dec, RA) axes.

Please checkout help(polseg_convert)
