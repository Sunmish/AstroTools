import numpy as np
import math

import astropy.units as u


c = 299792458.0

units = {"arcsec": u.arcsec, \
         "arcmin": u.arcmin, \
         "deg"   : u.degree, \
         "degree": u.degree, \
         "as"    : u.arcsec, \
         "am"    : u.arcmin, \
         "radian": u.radian, \
         "rad"   : u.radian
         }



def sBeam(freq, max_b, k=1.22, freq_unit=1.e6, return_unit="arcsec", \
	      precision=3):
	"""Calculate synthesized beam.

	Parameters
	----------
	freq        : float
				  Centre frequency of observation.
	max_b       : float
				  The maximum baseline in metres.
	k           : float, optional
				  Constant for resolution equation. Default is 1.22/'Airy disk'.
	freq_unit   : float, optional
				  Unit of frequency input. Default is MHz.
	return_unit : string or astropy.unit unit, optional
			      Units for return value. Default is arcsec.
	precision   : int, optional
				  Precision of final result.
 	"""

	s = k * c/(freq*freq_unit * max_b) * u.radian
	s = s.to(units[return_unit])

	pixel = s / 3.0

	print(">>> Synthesided beam : {0} {1}".format(\
		  round(s.value, precision), s.unit))
	print(">>> Pixel size       : {0} {1}".format(\
		  round(pixel.value, precision), pixel.unit))

	return s, pixel



def pBeam(freq, dish, k=1.22, freq_unit=1.e6, return_unit="arcsec", \
	      max_b=None, precision=3):
	"""Calculate primary beam size.\

	Parameters
	----------
	freq        : float
				  Centre frequency of observation.
	dish        : float
				  Dish diameter in metres.
	k           : float, optional
				  Constant for resolution equation. Default is 1.22/'Airy disk'.
	freq_unit   : float, optional
				  Unit of frequency input. Default is MHz.
	return_unit : string or astropy.unit unit, optional
			      Units for return value. Default is arcsec.
	max_b       : float, optional
				  The maximum baseline in metres. If provided, `pBeam` will 
				  calculate the synthesized beam and determine the suitable
				  number of pixels to cover the primary beam.
	precision   : int, optional
				  Precision of final result.
 	"""

	p = k * c/(freq*freq_unit * dish) * u.radian
	p = p.to(units[return_unit])

	print(">>> Primary beam     : {0} {1}".format(\
		  round(p.value, precision), p.unit))

	if max_b is not None:
		s, pixel = sBeam(freq, max_b, k, freq_unit, return_unit, precision)
		total_pixels = int(math.floor(p.value / pixel.value))

		print (">>> Total pixels     : {0} pixels".format(total_pixels))

		return p, total_pixels

	else:

		return p


def largest_scale(freq, min_b, k=1.22, freq_unit=1.e6, return_unit="arcsec", \
				  precision=3):
	"""Calculate the largest detectable angular scale given a minimum baseline.

	Parameters
	----------
	freq        : float
				  Centre frequency of observation.
	min_b       : float
				  The minimum baseline of the array in metres.
	k           : float, optional
				  Constant for resolution equation. Default is 1.22/'Airy disk'.
	freq_unit   : float, optional
				  Unit of frequency input. Default is MHz.
	return_unit : string or astropy.unit unit, optional
			      Units for return value. Default is arcsec.
	precision   : int, optional
				  Precision of final result.
 	"""

	l = k * c/(freq*freq_unit * min_b) * u.radian
	l = l.to(units[return_unit])

	print(">>> Largest detectable scale: {0} {1}".format(\
		  round(l.value, precision), l.unit))

	return l