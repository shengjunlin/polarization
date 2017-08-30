# polarization
Python functions related to ploarization.

In the python shell, import easy\_plot with the following code:  

    import sys
    sys.path.append('the path of the dir containing the script')
    from polarization import *

**polseg_convert()**

Generate a ds9 region file which contains polarization segments.

polseg_convert(I_map, polI_map, polPA_map, scale_10percent, sampling_rate, seg_color, output_reg)

Generate a ds9 region file which contains polarization segments.

Assume the fits files have 4 axes. (e.g. CASA simulation outputs)

Please checkout help(polseg_convert)
