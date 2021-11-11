#!/usr/bin/env python
"""Copy the relevant directories from validation_data_cfht, select only the
necessary detectors, zero out the imageHDUs in the calexps, and create any
necessary sylinks. Deletes `cfht/repo` before running: that path should only
contain the files that are handled by this script; refcats, export files, etc.
should live elsewhere.

No need to re-compress with fpack, as they are all already CompImageHDUs.
Requires that validation_data_cfht be setup, that this be run from
`testdata_jointcal/`.
"""
from astropy.io import fits
import os
import glob
import shutil
import sys

import pyarrow.parquet

import lsst.utils

# We only use these detectors for the jointcal tests.
detectors = (12, 13, 14, 21, 22, 23)

inpath = os.path.join(lsst.utils.getPackageDir("validation_data_cfht"), "repo")

# cleanup directories before copying new data
shutil.rmtree('cfht/repo', ignore_errors=True)

# copy new data
paths = glob.glob(f"{inpath}/singleFrame/*")
if len(paths) > 1:
    print(f"ERROR: Found {len(paths)} singleFrame output paths.")
    print("ERROR: Only run this script on a single-pipeline-run validation_data_cfht output repo.")
    sys.exit(-1)
run = paths[0].split('/')[-1]
shutil.copytree(f"{inpath}/singleFrame/{run}/calexp/",
                f"cfht/repo/singleFrame/{run}/calexp")
# TODO: do we need the config files?
# shutil.copytree(f"{inpath}/data/output/config/", "cfht/repo/config/")
shutil.copytree(f"{inpath}/singleFrame/{run}/sourceTable_visit/",
                f"cfht/repo/singleFrame/{run}/sourceTable_visit/")
shutil.copytree(f"{inpath}/singleFrame/{run}/visitSummary/",
                f"cfht/repo/singleFrame/{run}/visitSummary/")
shutil.copytree(f"{inpath}/skymaps", "cfht/repo/skymaps")
shutil.copytree(f"{inpath}/MegaPrime/calib", "cfht/repo/MegaPrime/calib")


def remove_detector_files(path):
    """Remove output files for detectors we don't need."""
    for file in glob.glob(path):
        ccd = int(file.split('ccd')[-1].split('_')[0])
        if ccd not in detectors:
            os.remove(file)


remove_detector_files(f"cfht/repo/singleFrame/{run}/calexp/*/r/r.MP9601/*/calexp*.fits")


def remove_detector_parquet(path):
    """Remove rows from parquet files for detectors we don't need."""
    for filename in glob.glob(path):
        data = pyarrow.parquet.ParquetFile(filename).read(use_pandas_metadata=True).to_pandas()
        match = data['detector'].isin(detectors)
        data[match].to_parquet(filename)


remove_detector_parquet(f"cfht/repo/singleFrame/{run}/sourceTable_visit/*/r/r.MP9601/*/*.parq")


# compress the calexps
files = glob.glob(f"cfht/repo/singleFrame/{run}/calexp/*/r/r.MP9601/*/calexp*.fits")
for f in files:
    data = fits.open(f)
    print('processing:', f)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU) or isinstance(hdu, fits.CompImageHDU):
            hdu.data[:] = 0
    data.writeto(f, overwrite=True)

# link in the refcats
os.symlink("../ref_cats", "cfht/repo/refcats")
