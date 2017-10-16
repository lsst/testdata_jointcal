#!/usr/bin/env python
"""Zero out the imageHDUs in the calexps, and then gzip the files."""
from __future__ import absolute_import, division, print_function

from astropy.io import fits
import glob
import subprocess

files = glob.glob('hsc/00*/HSC-?/corr/CORR-*.fits')
for f in files:
    data = fits.open(f)
    print('processing:', f)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU):
            hdu.data[:] = 0
    data.writeto(f, clobber=True)
    # fpack the file, using gzip for slightly better compression of all-zeros.
    # and '-qz 4' to force not dithering the zeros just in case.
    subprocess.call(['fpack', '-w', '-g2', '-qz', '4', f])
