#!/usr/bin/env bash

# This script will convert the hsc repo from gen2 to gen3.  This conversion requires the
# raws, so compress_jointcal_hsc_raw_test_data.py must have been run prior to this script.
# Until the gen3 repo format is fully stable, only these conversion scripts will be
# available in the repo, and not the generated gen3 repo.  After the gen3 format is
# stable enough, we can commit the output gen3 repo.

# Run the gen2 to gen3 conversion code
# This conversion is done in-place so the same path is used for the gen2 and gen3 repos.
butler convert hsc -i lsst.obs.subaru.HyperSuprimeCam --gen2root ./hsc

# Delete 99% of the skymap which is unused (and takes up almost 900 Mb)
sqlite3 hsc/gen3.sqlite3 "delete from patch_htm7_overlap where tract != 9697;" "delete from patch where tract != 9697;" "delete from tract_htm7_overlap where tract != 9697;" "vacuum;" ".exit"
