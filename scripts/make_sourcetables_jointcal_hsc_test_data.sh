#!/usr/bin/env sh

# Convert all source catalogs in the hsc/ repo to sourceTable_visit using
# default conversions.
# Should be run from the testdata_jointcal base path
# Note that --no-versions is specified because of the different version of the stack
# used to generate the source tables (w_2020_22) vs the original data (w_2020_14).
# This is necessary because we need obs_subaru from w_2020_18 or later
# (including DM-24379) in order to do the sourceTable conversions, and from
# w_2020_21 or later (including DM-25020) which includes the correct local background
# measurement.  The default Source.yaml file needs to be modified to add
# the sky_source column which did not exist for the w_2020_14 RC2 reprocessing.

# The doApplyExternalPhotoCalib and doApplyExternalSkyWcs options here refer to
# "external to the calexp src table", and by default will compute base_LocalPhotoCalib
# and base_LocalWcs from the calexp values unless otherwise overridden.  See
# DM-24379 for details.
writeSourceTable.py ./hsc --output ./hsc --id filter=HSC-R^HSC-I --config doApplyExternalPhotoCalib=True --config doApplyExternalSkyWcs=True --no-versions

transformSourceTable.py ./hsc --output ./hsc --id filter=HSC-R^HSC-I --no-versions --config functorFile=./scripts/policy/Source_hsc.yaml

consolidateSourceTable.py ./hsc --output ./hsc --id visit=34648^34670^34674^34690^34714^35892^36140^36192^36236^36260 --no-versions

# Finally, clean out SRC and sourceTable per-ccd parquet files that we do not need

rm hsc/0????/HSC-?/output/*-???.parq
