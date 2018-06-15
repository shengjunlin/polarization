# polarization
Python functions related to ploarization.

In the python shell, import polarization with the following code:  

    import sys
    sys.path.append('the path of the dir containing the script')
    from polarization import *

**polseg_convert()**

Sytnax:

    polseg_convert(I_map, polI_map, polPA_map, scale_10percent, sampling_interval=3, I_clip=0., polI_clip=0., seg_color='r', output_reg='output.reg', hist_plot=True)

Arguments:

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

Generate a ds9 region file which contains polarization segments.

Assume the fits files have 4 axes. (e.g. CASA simulation outputs)

Please checkout help(polseg_convert)
