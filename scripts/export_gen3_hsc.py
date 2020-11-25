#!/usr/bin/env python
# Create a gen3 exports.yaml file from an existing gen3 hsc repo.
# Requires that `convert_gen2_to_gen3_hsc.sh` has already been run.

import lsst.daf.butler as dafButler
from lsst.daf.butler import CollectionType
from lsst.obs.subaru import HyperSuprimeCam

butler = dafButler.Butler('hsc/butler.yaml')

collection = HyperSuprimeCam.makeCollectionName('testdata')

with butler.export(filename="hsc/exports.yaml") as export:
    # Raws
    rawCollection = HyperSuprimeCam.makeDefaultRawIngestRunName()
    export.saveDatasets(butler.registry.queryDatasets(collections=rawCollection,
                                                      datasetType='raw'))

    # Datasets
    export.saveDatasets(butler.registry.queryDatasets(collections=collection,
                                                      datasetType=...))

    # Calibrations
    for datasetTypeName in ('camera', 'transmission_optics', 'transmission_sensor',
                            'transmission_filter', 'transmission_atmosphere'):
        export.saveDatasets(butler.registry.queryDatasets(datasetTypeName, collections=...))

    for collection in butler.registry.queryCollections(..., collectionTypes={CollectionType.CALIBRATION}):
        export.saveCollection(collection)

    # Reference catalog
    export.saveDatasets(butler.registry.queryDatasets(collections='refcats',
                                                      datasetType='ps1_pv3_3pi_20170110'))
    export.saveDatasets(butler.registry.queryDatasets(collections='refcats',
                                                      datasetType='gaia_dr2_20200414'))

    # The skymap data
    export.saveDatasets(butler.registry.queryDatasets(collections='skymaps',
                                                      datasetType='deepCoadd_skyMap'))
    export.saveDataIds(
        butler.registry.queryDataIds(
            ["patch"],
            where="tract=9697 AND skymap='hsc_rings_v1'"
        ).expanded()
    )
