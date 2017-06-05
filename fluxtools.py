# TODO: 
# Add precision options. Currently values are returned at arbitrary precision.
# Move utility functions to separate file? 
#

import numpy
import math
import os
import sys

from datetime import datetime             # For timing purposes.

from astropy.io import fits               # For handling FITS files.
from astropy.wcs import WCS               # For computing LAS/positions.
from astropy.coordinates import SkyCoord  # For cross-referencing.
# SkyCoord or coordinates are returning import errors?
# This is needed for `measure_tree` - finding a source in a catalogue.
from astropy import units as u            

import logging
logging.basicConfig(format="%(levelname)s (%(module)s): %(message)s", \
                    level=logging.INFO)




__author__  = "Stefan Duchesne"
__version__ = "v1.0"
__date__    = "06-06-2017" 



units = {"arcsec": u.arcsec, \
         "arcmin": u.arcmin, \
         "deg"   : u.degree, \
         "degree": u.degree, \
         "as"    : u.arcsec, \
         "am"    : u.arcmin, \
         "radian": u.radian, \
         "rad"   : u.radian
         }



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




def order_boundaries(coords):
    """Order boundary coordinates based on nearest-neighbours.

    Only used for polygon annotations that are not too messy.
    """

    un_x = [coord[0] for coord in coords]
    un_y = [coord[1] for coord in coords]
    ord_coord = [(un_x.pop(0), un_y.pop(0))]

    for x, y in ord_coord:
        init  = True
        index = None
        for i in range(len(un_x)):
            d = numpy.sqrt((x-un_x[i])**2 + (y-un_y[i])**2)
            if init:
                diff  = d
                index = i
                init  = False
            elif d < diff:
                diff  = d
                index = i
            else:
                pass

        if index is None:
            ord_coord.append(ord_coord[0])
            break
        else:
            ord_coord.append((un_x.pop(index), un_y.pop(index)))

    
    return ord_coord



def kvis_polygon(region, r_index=0):
    """Get coordinate list for polygon (CLINES) in annotation file."""

    count = 0

    f = open(region, "r")
    lines = f.readlines()

    polygons = []
    for line in lines:
        if "clines" in line.lower():
            polygons.append(line)

    poly = polygons[r_index]  # If there are more than one cline entries.

    bits = poly.split(" ")
    coords = []
    for i in range(1, len(bits)-1, 2):
        coords.append((float(bits[i]), float(bits[i+1])))

    f.close()

    return coords
 


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
        raise ValueError(">>> NAXIS must be 2, 3, or 4.")

    if len(ra) == 1: ra, dec = ra[0], dec[0]

    return ra, dec



def world_to_pix(ra, dec, warray, naxis):
    """A wrapper for astropy's all_world2pix.

    This helps with issues where NAXIS > 2.

    Returns lists if the input pixel coordinates are lists/arrays.
    """

    if (not isinstance(ra, list)) or (not isinstance(ra, numpy.ndarray)):
        ra, dec = [ra], [dec]

    if naxis == 2:
        y, x = warray.all_world2pix(ra, dec, 0)
    elif naxis == 3:
        y, x = warray.all_world2pix(ra, dec, numpy.ones_like(ra), 0)
    elif naxis == 4:
         y, x = warray.all_world2pix(ra, dec, numpy.ones_like(ra), \
                                     numpy.ones_like(ra), 0)
    else:
        raise ValueError(">>> NAXIS must be 2, 3, or 4.")

    if len(x) == 1: x, y = x[0], y[0]

    return x, y



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


def rms_array(rms, farray):
    """Generate rms array or 

    """

    try:
        rms = float(rms)
        rarray = None
    except (TypeError, ValueError):
        if isinstance(rms, str): rarray = fits.getdata(rms)
        elif isinstance(rms, fits.HDUList): rarray = rms[0].data
        else: raise ValueError(">>> RMS must be specified as either a single " \
                               "value or as an array/filepath.")
    if rarray is None:
        rarray = numpy.full_like(farray, rms, dtype=numpy.double)

    if rarray.shape != farray.shape: 
        raise ValueError(">>> RMS array and image array must be the same size.")

    return rarray



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
            self.bounds = numpy.array([(m, n)])              # Boundary coordinates
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

                if i == m and j == n: boundary = True
                else: boundary = False

                for index in surrounding_indices:


                    if (index[0] < 0) or (index[1] < 0):
                        pass
                    else:
                        try:
                            if numpy.isnan(forest[index]) and \
                                (farray[index] >= self.gh*rarray[index]):
                                self.leaves += 1
                                self.fluxes = numpy.append(self.fluxes, \
                                    farray[index])
                                self.rluxes = numpy.append(self.rluxes, \
                                    rarray[index])
                                self.coords = numpy.append(self.coords, [index], \
                                    axis=0)
                                forest[index] = self.no
                                if farray[index] > self.bright:
                                    self.bright = farray[index]
                                    self.bcoord = [index]
                                indices.append(index)
                            elif numpy.isnan(forest[index]) and \
                                (farray[index] < self.gh*rarray[index]):
                                forest[index] = 0
                                farray[index] = numpy.nan
                                if not boundary:
                                    self.bounds = numpy.append(self.bounds, \
                                        [(i, j)], axis=0)
                                    boundary = True
                            elif forest[index] == 0:
                                if not boundary:
                                    self.bounds = numpy.append(self.bounds, \
                                        [(i, j)], axis=0)
                                    boundary = True
                                pass
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

    ann_opened = False

    if outfile is not None: start_time = datetime.now()
    # --------------------------------------------------------------------------
    if annfile is not None:
        if annfile.endswith("ann") or annfile.endswith("reg"): 
            ann = open(annfile, "w+")
            ann_opened = True
        else:
            logging.warn(">>> Annotation file must be for Kvis (.ann) or " \
                         " DS9 (.reg).")
            logging.warn(">>> No annotation file will be made.")
    # --------------------------------------------------------------------------


    if cutoff1 is None: cutoff1 = 3
    if cutoff2 is None: cutoff2 = cutoff1
    logging.info(">>> Detection threshold set to {0} sigma.".format(cutoff1))
    logging.info(">>> Growth threshold set to ...{0} sigma.".format(cutoff2))

    farray, warray, bpp, cd1, cd2, naxis = read_fits(fitsimage)
    rarray = rms_array(rms, farray)

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

    source_flux     = []
    source_dflux    = []
    source_peak     = []
    source_centroid = []
    source_avg_flux = []
    source_avg_rms  = []
    source_bcoord   = []
    source_area     = []
    source_npix     = []
    source_LAS      = []
    source          = []
    source_boundaries = []

    for tree in tree_leaves:

        source_flux.append(sum(tree_fluxes[tree]) / bpp)
        source_dflux.append((sum(tree_rms[tree]) / bpp) * numpy.sqrt(bpp / \
                            float(tree_leaves[tree])))
        source_avg_flux.append(sum(tree_fluxes[tree]) / tree_leaves[tree])
        source_avg_rms.append(sum(tree_rms[tree]) / tree_leaves[tree])
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

        # ----------------------------------------------------------------------
        if ann_opened:
            ord_coord = order_boundaries(tree_bounds[tree])
            if annfile.endswith("ann"):
                ann.write("COLOR RED\nCOORD W\n")
                ann.write("CLINES ")
            elif annfile.endswith("reg"):
                ann.write("global color=yellow\nfk5\n")
                ann.write("polygon ")
            for i in range(len(ord_coord)):
                r_dd, d_dd = pix_to_world(ord_coord[i][0], \
                                          ord_coord[i][1], warray, naxis)
                ann.write("{0} {1} ".format(r_dd, d_dd))
            ann.write("\n")
        # ----------------------------------------------------------------------


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
                    " # source{9}flux{9}dflux{9}avg_flux{9}avg_rms{9}area{9}npix{9}centroid_RA{9}centroid_DEC{9}bright_RA{9}bright_DEC{9}LAS\n" \
                    " #         Jy    Jy  Jy/Beam  deg^2   Jy/Beam   deg         deg          deg      deg     deg\n" \
                    " #   0     1     2      3      4     5       6           7            8        9       10   11\n" \
                    .format(fitsimage, outfile, cutoff1, cutoff2, min_pix, \
                    max_pix, len(source), datetime.now()-start_time, LAS, de))

            for i in range(len(source)):
                f.write("{0}{11}{1}{11}{2}{11}{3}{11}{4}{11}{5}{11}{6}{11}{7}" \
                        "{11}{8}{11}{9}{11}{10}{11}{12}\n".format(source[i], \
                        source_flux[i], source_dflux[i], source_avg_flux[i], \
                        source_avg_rms[i], source_area[i], source_npix[i], \
                        world_coords[0][i], world_coords[1][i], \
                        bright_coords[0][i], bright_coords[1][i], de, \
                        source_LAS[i]))

    # --------------------------------------------------------------------------
    if ann_opened:

        if annfile.endswith("reg"):
            for i in range(len(source)):
                ann.write(r"text {0} {1} {{0}}\\n".format(\
                          world_coords[0][i], world_coords[1][i], \
                          source[i]))

        elif annfile.endswith("ann"):
            for i in range(len(source)):
                ann.write("TEXT {0} {1} {2}\n".format( \
                          world_coords[0][i], world_coords[1][i], \
                          source[i]))

        ann.close()
    # --------------------------------------------------------------------------

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



def measure_tree(fitsimage, coords, rms, cutoff1=3, cutoff2=None, max_pix=500, \
                 min_pix=2, diagonals=True, LAS=True, annfile=None, outfile=None, \
                 outimage=None, verbose=False):
    """Calculate the fluxes of an individual.

    Parameters
    ----------
    fitsimage   : str or HDUList object
                Either the path to a FITS file or the HDUList object if already
                opened.
    coords      : tuple of floats
                Coordinates of source of interest. Should be in decimal degrees.
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

    # A forest in which our tree resides:
    source, source_flux, source_dflux, source_avg_flux, source_area, \
        source_npix, world_coords, bright_coords, source_LAS, source_peak = \
        measure_forest(fitsfile, rms, cutoff1, cutoff2, max_pix, min_pix, \
        diagonals, LAS, annfile, outfile, outimage, verbose)

    # Now to find the tree:
    c = SkyCoord(coords[0], coords[1], unit=(u.deg, u.deg))
    ww_catalogue = SkyCoord(world_coords, unit=(u.deg, u.deg))
    i = c.match_to_catalog_sky(ww_catalogue)[0]

    dist = angular_distance((coords[0], coords[1]), \
                            (world_coords[0][i], world_coords[1][i]))

    if verbose:
        if dist < 1.0:
            dist_print = dist * 60.0
            dist_unit = "arcmin"
            if dist_print < 1.0:
                dist_print *= 60.0
                dist_unit = "arcsec"
        else:
            dist_print = dist
            dist_unit = "deg"
        print(">>> The following parameters have been found:\n" \
              ".............. Source no. = {11}\n" \
              "............... Int. flux = {0} [Jy]\n" \
              "......... Error int. flux = {1} [Jy]\n" \
              "............... Peak flux = {2} [Jy/beam]\n" \
              "............... Avg. flux = {3} [Jy/beam]\n" \
              ".............. No. pixels = {4}\n" \
              "Flux weighted coordinates = ({5}, {6})\n" \
              "..................... LAS = {7} [deg]\n" \
              "............. Source area = {8} [deg^2]\n" \
              "....... Offset from input = {9} [{10}]".format(\
              source_flux[i], source_dflux[i], source_peak[i], \
              source_avg_flux[i], source_npix[i], world_coords[0][i], \
              world_coords[1][i], source_LAS[i], source_area[i], dist_print, \
              dist_unit, source[i]))


    return source[i], source_flux[i], source_dflux[i], source_peak[i], \
        source_avg_flux[i], source_npix[i], world_coords[i], source_LAS[i], \
        source_area[i], bright_coords[i]   



def measure_aperture(fitsimage, coords, radius, rms=None, sigma=3, LAS=True, \
                     verbose=True, radius_units="deg"):
    """Measure flux within a circular aperture.

    Parameters
    ----------
    fitsfimage   : string or astropy.io.fits.HDUList
                 If string this should be the filename and path to a FITS file.
    RA           : float
                 Central RA in decimal degrees.
    DEC          : float
                 Central DEC in decimal degrees.
    radius       : float
                 Aperture within which to calculate flux in degree.
     rms         : float or str
                 Either a single rms value for the image, or a rms image mirroring
                 the input FITS file. Minimum detection threshold is `rms` * `cutoff`
    cutoff       : int, optional
                 Multiple of the RMS required for detection. Default is 3.
    LAS          : bool, optional
                 If True an LAS is calculated.
    verbose      : bool, optional
                 If True results are printed to the terminal.
    radius_units : {'deg', 'arcmin', 'arcsec'}, optional
                 Specifies the unit for `radius`.

    """

    farray, warray, bpp, cd1, cd2, naxis = read_fits(fitsimage)
    rarray = rms_array(rms, farray)

    r = radius * units[return_unit]

    source_flux, source_rms, source_xpixel, source_ypixel, source_coords = \
     [], [], [], [], []
    source_peak = 0


    for i in range(len(farray[:, 0])):
        for j in range(len(farray[0, :])):

            if farray[i, j] >= cutoff*rarray[i, j]:

                c = pix_to_world(i, j, warray, naxis)
                diff = angular_distance(coords, (c[0], c[1]))

                if diff <= (radius.to(u.degree).value):

                    source_flux.append(farray[i, j])
                    source_rms.append(rarray[i, j])
                    source_xpixel.append(i)
                    source_ypixel.append(j)
                    source_coords.append(c)
                    if farray[i, j] > source_peak: source_peak = farray[i, j]

    int_flux = sum(source_flux) / bpp
    unc_flux = (sum(source_rms) / bpp) * numpy.sqrt(bpp / float(len(source_rms)))
    area     = len(source_flux) * abs(cd1*cd2)
    npix     = len(source_flux)

    if LAS:
        las = 0
        for i in range(len(source_coords)):
            for j in range(len(source_coords)):
                if angular_distance(source_coords[i], source_coords[j]) > las:
                    las = angular_distance(source_coords[i], source_coords[j])
    else:
        las = "NA"

    if verbose:
        print(">>> The following parameters have been found:\n" \
              "............... Int. flux = {0} [Jy]\n" \
              "......... Error int. flux = {1} [Jy]\n" \
              "............... Peak flux = {2} [Jy/beam]\n" \
              ".............. No. pixels = {3}\n" \
              "..................... LAS = {4} [deg]\n" \
              "............. Source area = {5} [deg^2]".format(\
              int_flux, unc_flux, source_peak, npix, las, area))


    return int_flux, unc_flux, source_peak, npix, las, area 




def measure_region(fitsimage, rms, region, r_index=0, sigma=3, annfile=None, \
                   annfile_color="yellow", verbose=True):
    """Measure integrated flux within polygon region. 

    Requires pyregion if using a ds9 region file. Inspired by `radioflux.py` by
    M.~J. Hardcastle: https://github.com/mhardcastle/radioflux
    
    This adapts a C++ implementation found here:
    http://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/


    Note that region vertices must be ordered.


    Parameters
    ----------
    fitsfimage   : string or astropy.io.fits.HDUList
                 If string this should be the filename and path to a FITS file.
    rms          : float or str
                 Either a single rms value for the image, or a rms image mirroring
                 the input FITS file. Minimum detection threshold is `rms` * `cutoff`
    region       : str, list, or numpy.ndarray
                 If a string, this should be a filepath to a ds9.reg file. If
                 a list or array, these should be ordered vertices of a polygon
                 in world coordinates of the image.
    r_index      : int, optional
                 Specifies the index of the polygon if using a ds9.reg or kvis.ann file. 
    sigma        : int, optional
                 Multiple of the RMS required for measurement. Default is 3.
    annfile      : str, optional
                 File to write annotations to. This creates a new file or over-
                 writes an existing one. The annotations here are just dots 
                 indicating each pixel that is measured.
    annfile_color: str, optional
                 Color of the annotations. Default is `yellow`.
    verbose      : bool, optional
                 If True, additional output is printed to the terminal.


    """




    def is_in(x, y, region_x, region_y, max_x):
        """Check if (x, y) is in region.

        (x, y) is considered in region if the followed is met:
        A line drawn horizontally to the right from (x, y) intersects 
        with an odd number of region edges. 
        """  
        

            
        def orientation(p1, p2, p3):
            """Get orientation of ordered triplet. 
            
            0 == collinear
            1 == clockwise
            2 == counter-clockwise
            """

            slope = (p2[1] - p1[1]) * (p3[0] - p2[0]) - \
                   (p3[1] - p2[1]) * (p2[0] - p1[0])

            if slope < 0:
                orient = 2
            elif slope > 0:
                orient = 1
            else:
                orient = 0


            return orient


        def on_segment(p1, p2, p3):
            """Checks if p3 lies on segment (p1, p2)."""

            if ((p3[0] <= max(p1[0], p2[0])) and (p3[0] >= min(p1[0], p2[0])) and\
                (p3[2] <= max(p1[1], p2[1])) and (p3[1] >= min(p1[1], p2[1]))):

                return True

            else:

                return False

        def intersects(p1, q1, p2, q2):
            """Determine if line segment (p1, q1) intersects with (p2, q2)

            Uses orientation conditions.
            """

            inter = False

            # General case:
            if (orientation(p1, q1, p2) != orientation(p1, q1, q2)) and \
               (orientation(p2, q2, p1) != orientation(p2, q2, q1)):

               inter = True

            # Special cases:
            if (orientation(p1, q1, p2) == 0) and (on_segment(p1, q1, p2)):
                inter = True
            if (orientation(p1, q1, q2) == 0) and (on_segment(p1, q1, q2)):
                inter = True
            if (orientation(p2, q2, p1) == 0) and (on_segment(p2, q2, p1)):
                inter = True
            if (orientation(p2, q2, q1) == 0) and (on_segment(p2, q2, q1)):
                inter = True


            return inter 

    
        coords_in_region = []

        for i in range(len(x)):

            p1, q1 = (x[i], y[i]), (x[i]+max_x, y[i])

            intersections = 0

            for j in range(len(region_x)):
                p2 = (region_x[j], region_y[j])
                if j == len(region_x)-1:
                    q2 = (region_x[0], region_y[0])
                else:
                    q2 = (region_x[j+1], region_y[j+1])

                if intersects(p1, q1, p2, q2):
                    intersections += 1

            if intersections%2 == 0:
                in_region = False
            else:
                in_region = True


            coords_in_region.append(in_region)

        return coords_in_region


    # 
    ra, dec = [], []

    if region.endswith(".reg"):

        import pyregion
        
        r = pyregion.open(region)[r_index].coord_list
        
        for i in range(0, len(r), 2):
            ra.append(r[i])
            dec.append(r[i+1])

    elif region.endswith(".ann"):
        
        poly_coords = kvis_polygon(region)
        ra  = [coord[0] for coord in poly_coords]
        dec = [coord[1] for coord in poly_coords]

    elif isinstance(region, str):
        raise TypeError(">>> `region` must be one of: DS9.reg, Kvis.ann, or a " \
                        " list coordinate tuples.")

    else:
        for vertex in region:
            ra.append(vertex[0])
            dec.append(vertex[1])


    farray, warray, bpp, cd1, cd2, naxis = read_fits(fitsimage)
    rarray = rms_array(rms, farray)
    max_x = len(farray[:, 0]+1)


    pix_r, pix_d = world_to_pix(ra, dec, warray, naxis)
    x_all, y_all, f_all, r_all = [], [], [], []
    for i in range(len(farray[:, 0])):
        for j in range(len(farray[0, :])):
            x_all.append(i)
            y_all.append(j)
            f_all.append(farray[i, j])
            r_all.append(rarray[i, j])

    cir = is_in(x_all, y_all, pix_r, pix_d, max_x)


    source_flux, source_rms, source_peak = [], [], 0
    final_ra, final_dec = [], []

    for i in range(len(cir)):
        if cir[i] and (f_all[i] >= sigma*r_all[i]):
            source_flux.append(f_all[i]), source_rms.append(r_all[i])
            ra_i, dec_i = pix_to_world(x_all[i], y_all[i], warray, naxis)
            final_ra.append(ra_i), final_dec.append(dec_i)
            if f_all[i] > source_peak:
                source_peak = f_all[i]

    if annfile is not None:
        if annfile.endswith(".ann"):
            with open(annfile, "w+") as g:
                g.write("colour {0}".format(annfile_color))
                for i in range(len(final_ra)):
                    g.write("DOT W {0} {1}\n".format(final_ra[i], final_dec[i]))

    int_flux = sum(source_flux) / bpp
    unc_flux = (sum(source_rms) / bpp) * numpy.sqrt(bpp / float(len(source_rms)))
    avg_rms  = sum(source_rms) / len(source_rms)
    area     = len(source_flux) * abs(cd1*cd2)
    npix     = len(source_flux)


    if verbose:
        print(">>> The following parameters have been found:\n" \
              "............... Int. flux = {0} [Jy]\n" \
              "......... Error int. flux = {1} [Jy]\n" \
              "............... Peak flux = {2} [Jy/beam]\n" \
              ".............. No. pixels = {3}\n" \
              # "..................... LAS = {4} [deg]\n" \
              "............. Source area = {4} [deg^2]\n" \
              "............. Average rms = {5} [Jy/beam]\n".format(\
              int_flux, unc_flux, source_peak, npix, area, avg_rms))


    return int_flux, unc_flux, source_peak, npix, area, avg_rms
