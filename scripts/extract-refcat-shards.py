#!/usr/bin/env python

"""Extract the relevant htm pixels for the Gaia DR2 and PS1 refcats on lsst-dev.

The output is written to each of the `ref_cats/` subdirectories under each
dataset (cfht, decam, hsc) in the currently setup testdata_jointcal dir.
"""

from astropy.coordinates import SkyCoord
import astropy.units as u
import os.path
import shutil

from esutil import htm

gaia_path = "/datasets/refcats/htm/v1/gaia_dr2_20200414"
ps1_path = "/datasets/refcats/htm/v1/ps1_pv3_3pi_20170110"

indexer = htm.HTM(depth=7)


def find_shards(center, radius):
    """Return the shards that intersect with the supplied circle."""
    return list(indexer.intersect(center.ra.to_value(u.degree), center.dec.to_value(u.degree),
                                  radius.to_value(u.degree), True))


def copy_files(instrument, shards, inpath):
    """Copy relevant files into `$TESTDATA_JOINTCAL_DIR/instrument/ref_cats`."""
    name = os.path.basename(inpath)
    outpath = os.path.join(os.path.expandvars("$TESTDATA_JOINTCAL_DIR"), instrument, "ref_cats", name)
    os.mkdir(outpath)
    shutil.copy(os.path.join(inpath, "config.py"), outpath)
    shutil.copy(os.path.join(inpath, "master_schema.fits"), outpath)
    for shard in shards:
        shutil.copy(os.path.join(inpath, f"{shard}.fits"), outpath)
    print(f"Copied {len(shards)} shard files from {name} into {outpath}")


def copy_all(instrument, center, radius):
    """Copy both refcats into the instrument repo."""
    shards = find_shards(center, radius)
    copy_files(instrument, shards, gaia_path)
    copy_files(instrument, shards, ps1_path)


center = SkyCoord(214.856821, +52.662694, frame='icrs', unit='deg')
radius = 0.3011315146269983 * u.degree
copy_all("cfht", center, radius)

center = SkyCoord(150.429282, +2.724456, frame='icrs', unit='deg')
radius = 0.6493234040345428 * u.degree
copy_all("decam", center, radius)

center = SkyCoord(320.876455, -0.309242, frame='icrs', unit='deg')
radius = 0.31055727581556886 * u.degree
copy_all("hsc", center, radius)
