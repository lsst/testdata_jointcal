#!/usr/bin/env python
"""Copy the relevant directories from RC/w_2020_14/DM-24359-sfm, zero out the
imageHDUs in the calexps, and then gzip the files.

No need to re-compress with fpack, as they are all already CompImageHDUs.
Requires that this be run on `lsst-dev`, and run from `testdata_jointcal/`, with
lsst_distrib setup via eups.
"""
import os
import shutil
import glob
import subprocess

from astropy.io import fits
import sqlite3

import lsst.daf.persistence

# This is the data identified by Eli Rykoff as suitable:
dataIds = (# Filter HSC-R
           {'visit': 34648, 'ccds': [51, 59, 67]},
           {'visit': 34690, 'ccds': [56, 64, 48]},
           {'visit': 34714, 'ccds': [19, 13, 18]},
           {'visit': 34674, 'ccds': [28, 21, 15]},
           {'visit': 34670, 'ccds': [92, 97, 93]},
           # Filter HSC-I
           {'visit': 36140, 'ccds': [51, 59, 67]},
           {'visit': 35892, 'ccds': [39, 31, 23]},
           {'visit': 36192, 'ccds': [20, 14, 8]},
           {'visit': 36260, 'ccds': [13, 7, 2]},
           {'visit': 36236, 'ccds': [93, 98, 87]})

inpath = "/datasets/hsc/repo/rerun/RC/w_2020_14/DM-24359-sfm"
butler = lsst.daf.persistence.Butler(inpath)

# cleanup directories before copying new data
shutil.rmtree('hsc/', ignore_errors=True)

os.mkdir('hsc')
# Write an old-school gen2 `_mapper` file; we can't use repositoryCfg.yaml
# because we won't have the raw parent repo on machines other than lsst-dev.
with open('hsc/_mapper', 'w') as file:
    file.write('lsst.obs.hsc.HscMapper\n')

# copy new data
shutil.copytree(f"{inpath}/config/", "hsc/config/")
shutil.copytree(f"{inpath}/schema/", "hsc/schema/")
shutil.copytree(f"{inpath}/deepCoadd/", "hsc/deepCoadd/")
# get the transmission curves from the calib repo
shutil.copytree(f"/datasets/hsc/calib/20200115/transmission/", "hsc/transmission/")

# create sqlite registry
conn = sqlite3.connect("hsc/registry.sqlite3")
cmd = "CREATE TABLE raw (id integer primary key autoincrement, taiObs text,expId text,pointing int,dataType text,visit int,dateObs text,frameId text,filter text,field text,pa double,expTime double,ccdTemp double,ccd int,proposal text,config text,autoguider int, unique(visit,ccd));"
conn.execute(cmd)
cmd = "CREATE TABLE raw_visit (visit int,field text,filter text,dateObs text,taiObs text, unique(visit));"
conn.execute(cmd)
cursor = conn.cursor()
cursor.execute('ATTACH DATABASE "/datasets/hsc/repo/registry.sqlite3" AS original')
for partial in dataIds:
    cmd = f"INSERT INTO raw_visit SELECT * FROM original.raw_visit WHERE visit={partial['visit']}"
    print('EXECUTING:', cmd)
    cursor.execute(cmd)
    for ccd in partial['ccds']:
        cmd = f"INSERT INTO raw SELECT * FROM original.raw WHERE visit={partial['visit']} AND ccd={ccd}"
        print('EXECUTING:', cmd)
        cursor.execute(cmd)
conn.commit()
conn.close()


def copy_data(name):
    """Copy the files for name (e.g. 'calexp', 'src') from inpath to `hsc/`"""
    for partial in dataIds:
        for ccd in partial['ccds']:
            dataId = {'visit': partial['visit'], 'ccd': ccd}
            path = butler.get(f'{name}_filename', dataId=dataId)[0]
            newpath = os.path.join('hsc', path.split('DM-24359-sfm/')[-1])
            if not os.path.exists(os.path.dirname(newpath)):
                os.makedirs(os.path.dirname(newpath))
            print(f'copying: {path} -> {newpath}')
            shutil.copy(path, newpath)


copy_data('src')
copy_data('calexp')


# compress the calexps
files = glob.glob('hsc/0129?/HSC-?/corr/CORR*.fits')
for file in files:
    data = fits.open(file)
    print('processing:', file)
    for hdu in data:
        if isinstance(hdu, fits.ImageHDU) or isinstance(hdu, fits.CompImageHDU):
            hdu.data[:] = 0
    data.writeto(file, overwrite=True)

# compress the source catalogs
files = glob.glob('hsc/0129?/HSC-?/output/SRC*.fits')
for file in files:
    subprocess.call(['gzip', file])
    print('gzipped catalog:', file)
