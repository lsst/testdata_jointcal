#!/usr/bin/env python
"""Zero out the imageHDUs in the calexps, and then gzip the files."""
from __future__ import absolute_import, division, print_function

from astropy.io import fits
import glob
import subprocess

files = glob.glob('decam/01768??/calexp/calexp*.fits')
for f in files:
    data = fits.open(f)
    print('processing:', f)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU):
            hdu.data[:] = 0
    data.writeto(f, clobber=True)
    # gzip the file
    subprocess.call(['gzip', f])
    # move it back to the original .fits filename to keep the butler happy.
    subprocess.call(['mv', f + '.gz', f])
