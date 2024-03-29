from lsst.daf.butler import Butler
from lsst.daf.butler import CollectionType
from lsst.obs.lsst import Latiss


butler = Butler("/repo/embargo")

collection = "u/erykoff/LATISS/2023-08A-07AB-05AB/testdata_subset"
instrument = "LATISS"
skymap = "latiss_v1"

visits = [
    2023051100320,
    2023051100357,
    2023051100390,
    2023051100406,
    2023051100448,
    2023051100454,
    2023051100278,
    2023051100473,
    2023051100263,
    2023051100509,
    2023051100304,
    2023051100431,
    2023051100547,
    2023051100379,
    2023051100495,
    2023051100489,
    2023051100401,
    2023051100280,
    2023051100303,
    2023051100508,
]

tracts = [5614, 5615]

with butler.export(
        filename="latiss_export_from_embargo_premod.yaml",
        directory="latiss/testdata",
        transfer="copy",
) as export:
    # visitSummary
    export.saveDatasets(
        butler.registry.queryDatasets(
            collections=collection,
            datasetType="visitSummary",
            where="visit in (visits)",
            bind={"visits": visits},
            instrument=instrument,
            findFirst=True,
        )
    )

    # isolated_star_cat
    export.saveDatasets(
        butler.registry.queryDatasets(
            collections=collection,
            datasetType="isolated_star_cat",
            where="tract in (tracts)",
            bind={"tracts": tracts},
            instrument=instrument,
            skymap=skymap,
            findFirst=True,
        )
    )

    # isolated_star_sources
    export.saveDatasets(
        butler.registry.queryDatasets(
            collections=collection,
            datasetType="isolated_star_sources",
            where="tract in (tracts)",
            bind={"tracts": tracts},
            instrument=instrument,
            skymap=skymap,
            findFirst=True,
        )
    )

    # Calibrations
    for datasetTypeName in (
            "camera",
            "transmission_filter",
            "transmission_optics",
            "transmission_sensor",
    ):
        export.saveDatasets(
            butler.registry.queryDatasets(
                datasetTypeName,
                collections=f"{instrument}/calib",
                instrument=instrument,
            )
        )

    # Reference catalog
    export.saveDatasets(
        butler.registry.queryDatasets(
            "atlas_refcat2_20220201",
            collections=["refcats"],
            where="tract in (tracts)",
            bind={"tracts": tracts},
            skymap=skymap,
        )
    )

    # Skymap data
    export.saveDatasets(
        butler.registry.queryDatasets(
            collections="skymaps",
            datasetType="skyMap",
            skymap=skymap,
        )
    )

    # And the calibration collection.
    collections = butler.registry.queryCollections(
        "LATISS/calib",
        flattenChains=True,
        includeChains=True,
    )
    for collection in collections:
        export.saveCollection(collection)

    collections = butler.registry.queryCollections(
        "LATISS/calib/unbounded",
        flattenChains=True,
        includeChains=True,
    )
    for collection in collections:
        export.saveCollection(collection)

# Now modify the runs; we need the full runs with timestamps here to make sure
# the testdata is labeled correctly.

preopsRuns = [
    "u/erykoff/LATISS/2023-08A-07AB-05AB/testdata_subset/20230824T164324Z",
    "LATISS/runs/AUXTEL_DRP_IMAGING_2023-08A-07AB-05AB/w_2023_33/PREOPS-3613/20230822T175428Z",
    "LATISS/runs/AUXTEL_DRP_IMAGING_2023-08A-07AB-05AB/w_2023_33/PREOPS-3613/20230821T153036Z"
]

testdataRun = "LATISS/testdata"

with open("latiss/testdata/latiss_export_from_embargo_premod.yaml", "r") as infile:
    with open("latiss/exports.yaml", "w") as outfile:
        for line in infile:
            for preopsRun in preopsRuns:
                if f"run: {preopsRun}" in line or f"name: {preopsRun}" in line:
                    line = line.replace(preopsRun, testdataRun)

            outfile.write(line)
