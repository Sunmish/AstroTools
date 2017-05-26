# AstroTools
A set of Python tools for mostly radio astronomy uses.



A note:

``fluxtools.py`` calculates integrated flux densities (along with associated errors) in the same manner as the source-finding software [Duchamp](http://www.atnf.csiro.au/people/Matthew.Whiting/Duchamp/) developed by Matthew Whiting. In fact, the original purpose of the ``measure_forest.py`` programme was to imitate Duchamp naively in its flux density measurements (not the source-finding part, in ``fluxtools.py`` source-finding is done in the most basic way possible) but allow the output of both flux-weighted coordinates for sources as well as brightest-pixel coordinates without having to run the programme twice. Further, I wanted a ``python`` implementation to use with other ``python`` processes I was running. 

The actual source-finding is more similar to the first part of [Aegean](https://github.com/PaulHancock/Aegean) by Paul Hancock in that it uses a simple flood-fill to detect sources and grow the detections. The ``measure_forest.py``script for command-line use also takes very generously from Paul Hancock's command-line implementation of [BANE](https://github.com/PaulHancock/Aegean/wiki/BANE) which is both an attempt to learn how to make ``python`` programmes easily usable from the command-line thanks to the ``argparse`` module but also to make it more easily to people who are not constantly working in an ``iPython`` console (or similar). Finally, ``measure_forest`` has support for the rms maps generated by BANE, though one can still supply a single rms value for the entire image if needed.
