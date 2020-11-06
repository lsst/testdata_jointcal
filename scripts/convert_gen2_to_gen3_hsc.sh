#!/usr/bin/env bash

set -e

# This script will convert the hsc repo from gen2 to gen3.  This conversion requires the
# raws, so compress_jointcal_hsc_raw_test_data.py must have been run prior to this script.
# Until the gen3 repo format is fully stable, only these conversion scripts will be
# available in the repo, and not the generated gen3 repo.  After the gen3 format is
# stable enough, we can commit the output gen3 repo.

# Run the gen2 to gen3 conversion code
# This conversion is done in-place so the same path is used for the gen2 and gen3 repos.

# Clean out old gen3 repo if it is there
if [ -f "hsc/butler.yaml" ]; then
    echo "Clearing out old butler.yaml"
    rm hsc/butler.yaml
fi

if [ -f "hsc/gen3.sqlite3" ]; then
    echo "Clearing out old gen3.sqlite3"
    rm hsc/gen3.sqlite3
fi

if [ -d "hsc/HSC/calib" ]; then
    echo "Clearing out old curated calibrations"
    rm -r hsc/HSC/calib
fi

# Do the basic conversion
butler convert hsc --gen2root ./hsc -C scripts/config/convertRepoHsc.py
# Delete the bfKernel and defects which we do not need
rm -r hsc/HSC/calib/unbounded/bfKernel
rm -r hsc/HSC/calib/curated/*T*/defects

# Delete 99% of the skymap which is unused (and takes up almost 1 Gb)
sqlite3 hsc/gen3.sqlite3 "delete from patch_skypix_overlap where tract != 9697;" "delete from patch where tract != 9697;" "delete from tract_skypix_overlap where tract != 9697;" "vacuum;" ".exit"
