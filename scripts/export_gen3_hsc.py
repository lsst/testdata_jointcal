#!/usr/bin/env python

import lsst.daf.butler as dafButler
from lsst.daf.butler import CollectionType
from lsst.obs.subaru import HyperSuprimeCam
import re

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
    for datasetTypeName in ('camera', 'transmission_optics',
                            'transmission_filter', 'transmission_atmosphere'):
        export.saveDatasets(butler.registry.queryDatasets(datasetTypeName, collections=...))

    for collection in butler.registry.queryCollections(..., collectionTypes={CollectionType.CALIBRATION}):
        export.saveCollection(collection)

    # Reference catalog
    export.saveDatasets(butler.registry.queryDatasets(collections='refcats',
                                                      datasetType='ps1_pv3_3pi_20170110'))

    # The skymap data
    export.saveDatasets(butler.registry.queryDatasets(collections='skymaps',
                                                      datasetType='deepCoadd_skyMap'))
    export.saveDataIds(
        butler.registry.queryDataIds(
            ["patch"],
            where="tract=9697"
        ).expanded()
    )
