#! /usr/bin/env python


import sys
import os
from argparse import ArgumentParser

from AstroTools import fluxtools


__author__ = "Stefan Duchesne"

if __name__ == "__main__":

    usage = "{0} [options] FilePath.fits".format(os.path.basename(__file__))

    parser = ArgumentParser(usage=usage, description="Calculate the flux " \
                            "densities of sources within a FITS file.")
    parser.add_argument("fitsimage", metavar="", help="")
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="", \
                       help="File to write out results to. This creates a new " \
                       "file or overwrites an exisiting one. Default is to not" \
                       " write out results.")
    parser.add_argument("-r", "--rms", dest="rms", metavar="", \
                       help="Either a single rms value for the input image, " \
                       "or an rms image mirroring the input FITS image file.")
    parser.add_argument("--cutoff1", dest="cutoff1", default=3, metavar="", \
                       type=float, help="The multiple of sigma required for " \
                       "detection of a source. Default = 3 sigma.")
    parser.add_argument("--cutoff2", dest="cutoff2", type=float, metavar="", \
                       help="The multiple of sigma required for growing a " \
                       "source detection. Default = cutoff1 * sigma.")
    parser.add_argument("--max_pix", dest="max_pix", default=500, metavar="", \
                       type=int, help="Maximum number of pixels in a detection. " \
                       "This is useful to set lower or higher depending on your " \
                       "image and whether time is of concern. Default = 500.")
    parser.add_argument("--min_pix", dest="min_pix", default=2, metavar="",\
                       type=int, help="Minimum number of pixels required " \
                       "for a detection. Default = 2.")
    parser.add_argument("--diagonals", dest="diagonals", default=True, metavar="",\
                       help="Specifies whether source detection concerns " \
                       "diagonally connected pixels or not. Default = True.")
    parser.add_argument("--LAS", dest="LAS", default=True, metavar="",\
                       help="Specifies whether or not to calculate the largest" \
                       " angular scale of the detected sources. This can be slow." \
                       " Default = True.")
    parser.add_argument("--annfile", dest="annfile", metavar="", \
                       help="File to write" \
                       " annotations to. This creates a new file or overwrites " \
                       "and existing one. Default is to not create such a file.")
    parser.add_argument("--outimage", dest="outimage", metavar="", \
                       help="This allows" \
                       " writing out a FITS file with pixels not part of detected" \
                       " sources set to NaN. Default is to not create such a file.")
    parser.add_argument("-v", "--verbose", dest="verbose", metavar="", \
                       default=False, help="Print additional lines of output " \
                       "while running. Default = False.")

    args = parser.parse_args()
    if args.fitsimage is None:
        parser.print_help()
        sys.exit()

    logging = fluxtools.logging
    logging.info(">>> Running {0} {1} {2}.".format(os.path.basename(__file__), \
                                               fluxtools.__version__, \
                                               fluxtools.__date__))

    if not os.path.exists(args.fitsimage):
        logging.error(">>> {0} not found.".format(args.fitsimage))
        logging.error(">>> Exiting.")
        sys.exit(1)
    if args.rms is None:
        logging.error(">>> RMS must be specified.")
        logging.error(">>> Exiting.")
        sys.exit(1)

    fluxtools.measure_forest(args.fitsimage, args.rms, args.cutoff1, \
                             args.cutoff2, args.max_pix, args.min_pix, \
                             args.diagonals, args.LAS, args.annfile, \
                             args.outfile, args.outimage, args.verbose)

