#!/usr/bin/env python
"""Copy in the relevant raw files from /datasets/hsc/repo, and zero out the
primaryHDUs and imageHDUs in the raws.

This script must be run on `lsst-dev`, and run from `testdata_jointcal/`, with
lsst_distrib setup via eups.  The result of this script was committed as part
of DM-24300.
"""

import os
import shutil
import subprocess

from astropy.io import fits

import lsst.daf.persistence as dafPersist
import lsst.utils


dataDir = lsst.utils.getPackageDir('testdata_jointcal')
repo = os.path.join(dataDir, 'hsc')
sourcerepo = '/datasets/hsc/repo'

butlerTest = dafPersist.Butler(repo)
butlerSource = dafPersist.Butler(sourcerepo)

subset = butlerTest.subset('src')

for dataRef in subset:
    # Confirm that this src catalog actually exists in the repo.
    if not butlerTest.datasetExists('src', dataId=dataRef.dataId):
        continue

    sourceRawUri = butlerSource.getUri('raw', dataId=dataRef.dataId)
    testRawUri = butlerTest.getUri('raw', dataId=dataRef.dataId, write=True)

    testPath = os.path.dirname(testRawUri)
    if not os.path.isdir(testPath):
        os.makedirs(testPath)
    if not os.path.isfile(testRawUri):
        shutil.copy(sourceRawUri, testRawUri)

    data = fits.open(testRawUri)
    print('processing: ', testRawUri)
    for hdu in data:
        if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)):
            hdu.data[:] = 0
    data.writeto(testRawUri, overwrite=True)
    subprocess.call(['fpack', testRawUri])
    os.remove(testRawUri)
    os.rename(testRawUri + '.fz', testRawUri)
