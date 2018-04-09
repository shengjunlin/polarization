# polarization
Python functions related to ploarization.

In the python shell, import polarization with the following code:  

    import sys
    sys.path.append('the path of the dir containing the script')
    from polarization import *

**polseg_convert()**

Sytnax:

    polseg_convert(I_map, polI_map, polPA_map, scale_10percent, sampling_interval, seg_color, output_reg)

Arguments:

I_map             [str]: The fits filename of Stokes I.    
polI_map          [str]: The fits filename of Polarized intensity.    
polPA_map         [str]: The fits filename of PA[deg] of polarization segments.    
scale_10percent [float]: Length[arcsec] of 10% polarization segments.    
sampling_interval [int]: (1/sampling rate of segments)[pixel].    
seg_color         [str]: Color of segments.    
output_reg        [str]: Output region filename.

Generate a ds9 region file which contains polarization segments.

Assume the fits files have 4 axes. (e.g. CASA simulation outputs)

Please checkout help(polseg_convert)
