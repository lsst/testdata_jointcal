#!/usr/bin/env python
"""Copy the relevant directories from validation_data_decam, zero out the
imageHDUs in the calexps, and then gzip the files.

No need to re-compress with fpack, as they are all already CompImageHDUs (the
stack does image compression when saving the Exposures).
Requires that `validation_data_decam` be setup, and that this be run from
`testdata_jointcal/`.
"""
from astropy.io import fits
import os
import shutil
import glob
import subprocess

import lsst.utils

inpath = lsst.utils.getPackageDir("validation_data_decam")

# cleanup directories before copying new data
shutil.rmtree("decam/0176837", ignore_errors=True)
shutil.rmtree("decam/0176846", ignore_errors=True)
shutil.rmtree("decam/config", ignore_errors=True)
shutil.rmtree("decam/schema", ignore_errors=True)

# copy new data
shutil.copytree(f"{inpath}/data/output/config/", "decam/config/")
shutil.copytree(f"{inpath}/data/output/schema/", "decam/schema/")
os.makedirs("decam/0176837")
os.makedirs("decam/0176846")
shutil.copytree(f"{inpath}/data/output/0176837/calexp/", "decam/0176837/calexp/")
shutil.copytree(f"{inpath}/data/output/0176837/src/", "decam/0176837/src/")
shutil.copytree(f"{inpath}/data/output/0176846/calexp/", "decam/0176846/calexp/")
shutil.copytree(f"{inpath}/data/output/0176846/src/", "decam/0176846/src/")

# Cleanup unneeded files
os.remove(f"decam/schema/icSrc.fits")

# NOTE: CCD 13 has nonsense astrometry in both visits
files = glob.glob("decam/01768??/*/*_13.fits")
for file in files:
    os.remove(file)

# compress the calexps
files = glob.glob('decam/01768??/calexp/calexp-*.fits')
for f in files:
    data = fits.open(f)
    print('processing:', f)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU) or isinstance(hdu, fits.CompImageHDU):
            hdu.data[:] = 0
    data.writeto(f, overwrite=True)

# compress the source catalogs
files = glob.glob('decam/01768??/src/*.fits')
for file in files:
    subprocess.call(['gzip', file])
    print('gzipped catalog:', file)
