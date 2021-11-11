#!/usr/bin/env python
# Create a gen3 exports.yaml file from an existing gen3 cfht repo.
# Requires that `compress_jointcal_cfht_test_data.py` has already been run.

import lsst.daf.butler as dafButler
import os.path

butler = dafButler.Butler(os.path.expandvars('$VALIDATION_DATA_CFHT_DIR/repo/butler.yaml'))

with butler.export(filename="cfht/exports.yaml") as export:
    # Datasets: we only want a few detectors.
    where = "instrument='MegaPrime' and detector in (12,13,14,21,22,23)"
    for datasetType in ('calexp', 'sourceTable_visit', 'visitSummary'):
        export.saveDatasets(butler.registry.queryDatasets(collections='singleFrame',
                                                          datasetType=datasetType,
                                                          where=where))

    # Calibrations
    export.saveDatasets(butler.registry.queryDatasets('camera', collections=...))

    # Reference catalogs are handled by ingest_files using the .ecsv files in cfht/

    # The skymap data: this uses a discrete skymap, so only one tract (0).
    export.saveDatasets(butler.registry.queryDatasets(collections='skymaps',
                                                      datasetType='skyMap'))
    export.saveDataIds(
        butler.registry.queryDataIds(
            ["patch"],
            where="tract=0 AND skymap='discrete'"
        ).expanded()
    )

    # The CHAINED collection we will use in tests.
    export.saveCollection('singleFrame')
    # Empty collections that are part of the above CHAINED collection.
    export.saveCollection('MegaPrime/raw/all')
    export.saveCollection('MegaPrime/calib/unbounded')
