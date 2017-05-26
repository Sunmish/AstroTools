# TODO: 
# Add precision options. Currently values are returned at arbitrary precision.
#

import numpy
import math
import os
import sys

from datetime import datetime             # For timing purposes.

from astropy.io import fits               # For handling FITS files.
from astropy.wcs import WCS               # For computing LAS/positions.
# from astropy.coordinates import SkyCoord  # For cross-referencing.
# SkyCoord or coordinates are returning import errors?
# This is needed for `measure_tree` - finding a source in a catalogue.
from astropy import units as u            

import logging
logging.basicConfig(format="%(levelname)s (%(module)s): %(message)s", \
                    level=logging.INFO)




__author__  = "Stefan Duchesne"
__version__ = "v1.0"
__date__    = "26-05-2017" 





def angular_distance(coords1, coords2):
    """Get the angular distance between a set of RA, DEC coordinates in [dd]."""

    cos_angle = math.sin(math.radians(coords1[1])) * \
                math.sin(math.radians(coords2[1])) + \
                math.cos(math.radians(coords1[1])) * \
                math.cos(math.radians(coords2[1])) * \
                math.cos(math.radians(coords1[0]) - math.radians(coords2[0]))

    try:
        gamma = math.degrees(math.acos(cos_angle))
    except ValueError:
        gamma = math.degrees(math.acos(min(1, max(cos_angle, -1))))

    return gamma



def pix_to_world(x, y, warray, naxis):
    """A wrapper for astropy's all_pix2world.

    This helps with issues where NAXIS > 2. 

    Returns lists if the input pixel coordinates are lists/arrays.
    """

    if (not isinstance(x, list)) or (not isinstance(x, numpy.ndarray)):
        x, y = [x], [y] 

    if naxis == 2:
        ra, dec = warray.all_pix2world(y, x, 0)
    elif naxis == 3:
        ra, dec = warray.all_pix2world(y, x, numpy.ones_like(x), 0)
    elif naxis == 4:
        ra, dec = warray.all_pix2world(y, x, numpy.ones_like(x), \
                                       numpy.ones_like(x), 0)
    else:
        raise ValueError(">>> NAXIS size must be 2, 3, or 4.")

    if len(ra) == 1: ra, dec = ra[0], dec[0]

    return ra, dec




def read_fits(fitsimage):
    """Read a FITS file with astropy.io.fits and return relevant data.

    Because the FITS standard is in fact not completely standard there are
    a few little things that need to be checked to look for the necessary
    header data.
    """

    opened = False
    if isinstance(fitsimage, str):
        # if not fitsimage.endswith(".fits"): fitsimage += ".fits"
        hdulist = fits.open(fitsimage)
        opened = True
    elif isinstance(fitsimage, fits.HDUList):
        hdulist = fitsimage
    else:
        raise IOError(">>> Specified FITS image must an astropy.io.fits.HDUList" \
                      "object or a filepath.")

    harray = hdulist[0].header
    farray = hdulist[0].data
    warray = WCS(harray)

    if opened: hdulist.close()

    naxis = harray["NAXIS"]

    # Read some things from the FITS header and check some things:

    try: farray.shape = (harray["NAXIS2"], harray["NAXIS1"])
    except ValueError:  raise ValueError(">>> FITS files must be flat.")

    # Try to find pixel sizes:
    try: 
        cd1, cd2 = harray["CDELT1"], harray["CDELT2"]
    except KeyError:
        try:
            cd1, cd2 = harray["CD1_1"], harray["CD2_2"]
        except KeyError:
            raise KeyError(">>> Cannot find key for pixel sizes.")

    # Beam parameters:
    semi_a, semi_b, beam_area, beam_not_found = None, None, None, True
    if ("BMAJ" in harray.keys()) and ("BMIN" in harray.keys()):
        semi_a = 0.5 * harray["BMAJ"]  # Normal keys.
        semi_b = 0.5 * harray["BMIN"]
        beam_not_found = False
    elif ("CLEANBMJ" in harray.keys()) and ("CLEANBMN" in harray.keys()):
        semi_a = 0.5 * harray["CLEANBMJ"]  # Abnormal keys, e.g. VLSSr.
        semi_b = 0.5 * harray["CLEANBMN"]
        beam_not_found = False
    else:
        try:  # Check for NVSS:
            for i in range(len(harray["HISTORY"])):
                if "'NRAO VLA SKY SURVEY" in harra["HISTORY"][i]:
                    semi_a = 0.5 * 1.2500e-2
                    semi_b = 0.5 * 1.2500e-2
                    beam_not_found = False
            if beam_not-found: raise KeyError
        except KeyError:
            try:  # AIPS has a tell-tale mark:
                for i in range(len(harray["HISTORY"])):
                    if "BMAJ" in harray["HISTORY"][i]:
                        l = []
                        for t in harray["HISTORY"][i].split():
                            try: l.append(float(t))
                            except ValueError: pass
                        semi_a = 0.5 * l[0]
                        semi_b = 0.5 * l[1]
                        beam_not_found = False
                if beam_not_found: raise KeyError
            except (KeyError, IndexError):
                try:  # Check for SUMSS, which does NOT have the beam information.
                      # We will use the declination-dependent beam size based on
                      # the reference pixel of the FITS file. 
                    for i in range(len(harray["HISTORY"])):
                        if "Sydney University Molonglo Sky Survey (SUMSS)" in \
                            harray["HISTORY"][i]:
                            decl = math.radians(abs(harray["CRVAL1"]))
                            beam_area = numpy.pi * (45.0**2 / 3600.0**2) * \
                                        (1.0 / math.sin(decl))
                            beam_not_found = False
                    if beam_not_found: raise KeyError
                except KeyError:
                    raise KeyError(">>> Beam parameters not found.")

    if beam_area is None:
        beam_area = numpy.pi * semi_a * semi_b
    bpp = beam_area / (abs(cd1*cd2) * numpy.log(2.0))


    return farray, warray, bpp, cd1, cd2, naxis



class Tree():
    """Grow a tree with branches and leaves.

    Grow a source from neighbouring pixels.
    """


    def __init__(self, tree_number, threshold):
        """ 'The seed of a tree of giants.' """

        self.no = tree_number   # The source ID.
        self.th = threshold[0]  # Threshold above rms needed for detection.
        self.gh = threshold[1]  # Threshold above rms needed for growing the detection.


    def seedling(self, m, n, farray, warray, rarray, forest, diagonals):
        """Start from a seedling and grow a tree."""

        if farray[m, n] >= self.th*rarray[m, n]:       # Detection!
            forest[m, n] = self.no                     # Set ref to ID
            self.leaves = 1                            # Count pixels
            self.fluxes = numpy.array([farray[m, n]])  # Add flux
            self.coords = numpy.array([(m, n)])        # Add pixel coordinates
            self.bounds = numpy.array([(m, n)])        # Boundary coordinates
            self.bright = farray[m, n]                 # Brightest pixel
            self.bcoord = [(m, n)]                     # BP coordinates
            self.rluxes = numpy.array([rarray[m, n]])  # Add rmses

            indices = [(m, n)]  # Indices to check. This is added to.

            for i, j in indices:
                # Surrounding pixels:
                if diagonals:
                    surrounding_indices = [(i-1, j-1), (i-1, j), (i, j-1), \
                        (i+1, j-1), (i-1, j+1), (i+1, j), (i, j+1), (i+1, j+1)]
                else:
                    surrounding_indices = [(i-1, j), (i, j-1), (i+1, j), \
                        (i, j+1)]

                boundary = False

                for index in surrounding_indices:
                    if (index[0] < 0) or (index[1] < 0):
                        pass
                    else:
                        try:
                            if (numpy.isnan(forest[index])) and \
                                (farray[index] >= self.gh*rarray[index]):
                                self.leaves += 1
                                self.fluxes = numpy.append(self.fluxes, \
                                    farray[index])
                                self.coords = numpy.append(self.coords, [index], \
                                    axis=0)
                                forest[index] = self.no
                                if farray[index] > self.bright:
                                    self.bright = farray[index]
                                    self.bcoord = [index]
                                indices.append(index)
                            elif numpy.isnan(forest[index]):
                                forest[index] = 0
                                farray[index] = numpy.nan
                                if not boundary:
                                    self.bounds = numpy.append(self.bounds, \
                                        [index], axis=0)
                                    boundary = True
                            else:
                                if not boundary:
                                    self.bounds = numpy.append(self.bounds, \
                                        [index], axis=0)
                                    boundary = True
                        except IndexError:
                            pass

            tree_number = self.no + 1

        else:
            tree_number = self.no
            forest[m, n] = 0
            farray[m, n] = numpy.nan


        return farray, forest, tree_number



def populate_forest(farray, warray, rarray, threshold, max_pix, min_pix, \
                    diagonals):
    """Grow trees in a forest; find sources in a field."""

    # An empty forest:
    forest = numpy.empty((len(farray[:, 0]), len(farray[0, :])))
    forest[:] = numpy.nan

    tree_number = 1                 # The current source ID,
    tree_leaves = {tree_number: 0}  # its pixels,
    tree_fluxes = {tree_number: 0}  # its flux values,
    tree_coords = {tree_number: 0}  # its pixel coordinates,
    tree_bounds = {tree_number: 0}  # its boundary coordinates,
    tree_bright = {tree_number: 0}  # Source brightest pixel coordinates.
    tree_height = {tree_number: 0}  
    tree_rms    = {tree_number: 0}  # Sum of rms.

    for n in range(0, len(farray[0, :])):
        for m in range(0, len(farray[:, 0])):

            # If forest[m, n] is not NaN, it has already been checked.
            if numpy.isnan(forest[m, n]):
                t = Tree(tree_number, threshold)

                farray, forest, tree_number = t.seedling(m, n, farray, \
                    warray, rarray, forest, diagonals)

                try:
                    if (min_pix <= t.leaves <= max_pix):
                        tree_leaves[tree_number-1] = t.leaves
                        tree_fluxes[tree_number-1] = t.fluxes
                        tree_coords[tree_number-1] = t.coords
                        tree_bounds[tree_number-1] = t.bounds
                        tree_bright[tree_number-1] = t.bcoord
                        tree_height[tree_number-1] = t.bright
                        tree_rms[tree_number-1]    = t.rluxes

                    else:
                        pass
                except AttributeError:
                    pass


    return farray, forest, tree_leaves, tree_fluxes, tree_coords, tree_bounds, \
        tree_bright, tree_height, tree_rms



def measure_forest(fitsimage, rms=None, cutoff1=None, \
                   cutoff2=None, max_pix=500, min_pix=2, diagonals=True, \
                   LAS=True, annfile=None, outfile=None, outimage=None, \
                   verbose=True):
    """Calculate the fluxes of individual sources within a FITS file.

    Parameters
    ----------
    fitsimage   : str or HDUList object
                Either the path to a FITS file or the HDUList object if already
                opened.
    rms         : float or str
                Either a single rms value for the image, or a rms image mirroring
                the input FITS file. Minimum detection threshold is `rms` * `cutoff`
                for each pixel.
    cutoff1     : int, optional
                The multiple of `rms` needed for a detection. Default is 3sigma.
    cutoff2     : int, optional
                The multiple of `rms` needed for growing sources. Default is `cutoff1`.
    max_pix     : int, optional
                Maximum number of pixels in a detection. Useful only to save
                time ignoring large sources (e.g., Fornax A) as LAS calculations
                are ridiculously slow.
    min_pix     : int, optional
                Minimum number of pixels for a detection.
    diagonals   : bool, optional
                Specifies whether source detection counts pixels only connected
                diagonally. Select True for diagonal detection.
    LAS         : bool, optional
                Select True is wanting to calculate the largest angular size of
                each source. THIS IS VERY SLOW.
    annfile     : str, optional
                File to write annotations to. This creates a new file or over-
                writes an existing one.
    outfile     : str, optional
                File to write out results to. This creates a new file or over- 
                writes an existing one.
    outimage    : str, optional
                This allows writing a FITS file with the same header/size as 
                the input image but with pixels below the detection/growth 
                thresholds blanked. An exisiting file with this name will 
                be deleted.
    verbose     : bool, optional
                Select True if wanting to print output to terminal.

    """

    if outfile is not None: start_time = datetime.now()


    print rms
    try:
        rms = float(rms)
        rarray = None
    except TypeError:
        if isinstance(rms, str): rarray = fits.getdata(rms)
        elif isinstance(rms, fits.HDUList): rarray = rms[0].data
        else: raise ValueError(">>> RMS must be specified as either a single " \
                               "value or as an array/filepath.")


    if cutoff1 is None: cutoff1 = 3
    if cutoff2 is None: cutoff2 = cutoff1
    logging.info(">>> Detection threshold set to {0} sigma.".format(cutoff1))
    logging.info(">>> Growth threshold set to ...{0} sigma.".format(cutoff2))

    farray, warray, bpp, cd1, cd2, naxis = read_fits(fitsimage)

    if rarray is None:
        rarray = numpy.full_like(farray, rms, dtype=numpy.double)

    if rarray.shape != farray.shape: 
        raise ValueError(">>> RMS array and image array must be the same size.")

    farray, forest, tree_leaves, tree_fluxes, tree_coords, tree_bounds,\
        tree_bright, tree_height, tree_rms = populate_forest(farray=farray, \
        warray=warray, rarray=rarray, threshold=(cutoff1, cutoff2), \
        max_pix=max_pix, min_pix=min_pix, diagonals=diagonals)

    if (len(tree_leaves) == 1) and (tree_leaves[1] == 0):
        raise RuntimeError(">>> No sources detected for {0} with threshold = "\
                           "{1} sigma.".format(fitsfile, cutoff1))

    zero_flag = False
    if tree_leaves[1] == 0:
        zero_flag = True
        del tree_leaves[1]
        del tree_fluxes[1]
        del tree_coords[1]
        del tree_bounds[1]
        del tree_bright[1]
        del tree_height[1]

    source_flux = []
    source_dflux = []
    source_peak = []
    source_centroid = []
    source_avg_flux = []
    source_bcoord = []
    source_area = []
    source_npix = []
    source_LAS = []
    source = []

    for tree in tree_leaves:

        source_flux.append(sum(tree_fluxes[tree]) / bpp)
        source_dflux.append((sum(tree_rms[tree]) / bpp) * numpy.sqrt(bpp / \
                            float(tree_leaves[tree])))
        source_avg_flux.append(sum(tree_fluxes[tree]) / tree_leaves[tree])
        source_area.append(abs(cd1*cd2) * tree_leaves[tree])
        source_npix.append(tree_leaves[tree])
        source_peak.append(tree_height[tree])
        source_bcoord.append((tree_bright[tree][0][1], tree_bright[tree][0][0]))
        fx, fy = [], []
        fz = sum(tree_fluxes[tree])
        for i in range(len(tree_fluxes[tree])):
            fx.append(tree_coords[tree][i][0] * tree_fluxes[tree][i])
            fy.append(tree_coords[tree][i][1] * tree_fluxes[tree][i])
        source_centroid.append((sum(fx) / fz, sum(fy) / fz))

        # The largest angular scale of the source. 
        # This is simply calculated as the largest angular separation between
        # any two pixels that comprise the source.
        if LAS:
            length = 0
            for pc1 in range(len(tree_bounds[tree])):
                ra1, dec1 = pix_to_world(tree_bounds[tree][pc1][1], \
                                         tree_bounds[tree][pc1][0], \
                                         warray, naxis)
                for pc2 in range(pc1+1, len(tree_bounds[tree])):
                    ra2, dec2 = pix_to_world(tree_bounds[tree][pc2][1], \
                                             tree_bounds[tree][pc2][0], \
                                             warray, naxis)

                    if (ra1 == ra2) and (dec1 == dec2): 
                            diff = 0.0
                    else:
                        diff = angular_distance((ra1, dec1), (ra2, dec2))

                    if diff > length:
                        length = diff
            source_LAS.append(length)

        else:
            source_LAS.append("NA")

        if zero_flag: source.append(tree-1)  
        else: source.append(tree)

        wx, wy, bx, by = [], [], [], []
        for i in range(len(source)):
            wx.append(source_centroid[i][0])
            wy.append(source_centroid[i][1])
            bx.append(source_bcoord[i][0])
            by.append(source_bcoord[i][1])

        world_coords = pix_to_world(wx, wy, warray, naxis)
        bright_coords = pix_to_world(bx, by, warray, naxis)

    
    if outfile is not None:

        if outfile[-3:] == "csv": de = ","
        else: de = " " 

        with open(outfile, "w+") as f:

            f.write(" # Sources and their parameters found by `measure_forest`:\n"\
                    " #\n" \
                    " # ............ Input FITS = {0}\n" \
                    " # ................ Output = {1}\n" \
                    " # .... Detection treshold = {2}\n" \
                    " # ...... Growth threshold = {3}\n" \
                    " # ........ Minimum pixels = {4}\n" \
                    " # ........ Maximum pixels = {5}\n" \
                    " # Total number of sources = {6}\n" \
                    " # ............ Time taken = {7}\n" \
                    " # ................... LAS = {8}\n" \
                    " # \n" \
                    " # source{9}flux{9}dflux{9}avg_flux{9}area{9}npix{9}centroid_RA{9}centroid_DEC{9}bright_RA{9}bright_DEC{9}LAS\n" \
                    " #         Jy    Jy  Jy/Beam  deg^2         deg         deg          deg      deg     deg\n" \
                    " #   0     1     2      3      4     5       6           7            8        9       10\n" \
                    .format(fitsimage, outfile, cutoff1, cutoff2, min_pix, \
                    max_pix, len(source), datetime.now()-start_time, LAS, de))

            for i in range(len(source)):
                f.write("{0}{11}{1}{11}{2}{11}{3}{11}{4}{11}{5}{11}{6}{11}{7}" \
                        "{11}{8}{11}{9}{11}{10}\n".format(source[i], \
                        source_flux[i], source_dflux[i], source_avg_flux[i], \
                        source_area[i], source_npix[i], world_coords[0][i], \
                        world_coords[1][i], bright_coords[0][i], \
                        bright_coords[1][i], source_LAS[i], de))


    if annfile is not None:

        if annfile[-3:] == "reg":
            with open(annfile, "w+") as reg:
                reg.write("global color=yellow\nfk5\n")

                for i in range(len(source)):
                    reg.write(r"text {0} {1} {{0}}\\n".format(\
                              world_coords[0][i], world_coords[1][i], \
                              source[i]))

        elif annfile[-3:] == "ann":
            with open(annfile, "w+") as ann:
                ann.write("COLOR YELLOW\nCOORD W\n")

                for i in range(len(source)):
                    ann.write("TEXT {0} {1} {2}\n".format( \
                              world_coords[0][i], world_coords[1][i], \
                              source[i]))
        else:
            raise ValueError(">>> Annotations must be for 'ds9' or 'kvis' " \
                             "with extension '.reg' or '.ann' respectively.")


    if outimage is not None:

        if isinstance(fitsimage, str):
            hdulist = fits.open(fitsimage)
        elif isinstance(fitsimage, fits.HDUList):
            hdulist = fitsimage

        hdulist[0].data = farray
        if os.path.exists(outimage):
            logging.warn(">>> Deleting old {0}...".format(outimage))
            os.remove(outimage)
        hdulist[0].header["HISTORY"] = "Output image from `fluxtools.py`."
        hdulist.writeto(outimage)
        hdulist.close()



    return source, source_flux, source_dflux, source_avg_flux, source_area, \
        source_npix, world_coords, bright_coords, source_LAS, source_peak



# def measure_tree()
