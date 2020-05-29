#!/usr/bin/env python
"""Copy the relevant directories from validation_data_cfht, zero out the
imageHDUs in the calexps, and then gzip the files.

No need to re-compress with fpack, as they are all already CompImageHDUs.
Requires that validation_data_cfht be setup, that this be run from
`testdata_jointcal/`.
"""
from astropy.io import fits
import os
import shutil
import glob
import subprocess

import lsst.utils

inpath = lsst.utils.getPackageDir("validation_data_cfht")

# cleanup directories before copying new data
shutil.rmtree('cfht/calexp', ignore_errors=True)
shutil.rmtree('cfht/config', ignore_errors=True)
shutil.rmtree('cfht/src', ignore_errors=True)
shutil.rmtree('cfht/metadata', ignore_errors=True)
shutil.rmtree('cfht/schema', ignore_errors=True)

# copy new data
shutil.copytree(f"{inpath}/data/output/calexp/", "cfht/calexp/")
shutil.copytree(f"{inpath}/data/output/config/", "cfht/config/")
shutil.copytree(f"{inpath}/data/output/src/", "cfht/src/")
shutil.copytree(f"{inpath}/data/output/schema/", "cfht/schema/")

# Cleanup unneeded files
os.remove(f"cfht/schema/icSrc.fits")
for file in glob.glob("cfht/calexp/06AL01/D3/2006-0*/r/bkgd*.fits"):
    os.remove(file)

# compress the calexps
files = glob.glob('cfht/calexp/06AL01/D3/2006-0*/r/calexp*.fits')
for f in files:
    data = fits.open(f)
    print('processing:', f)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU) or isinstance(hdu, fits.CompImageHDU):
            hdu.data[:] = 0
    data.writeto(f, overwrite=True)

# compress the source catalogs
files = glob.glob('cfht/src/06AL01/D3/2006-0*/r/SRC*.fits')
for file in files:
    subprocess.call(['gzip', file])
    print('gzipped catalog:', file)
