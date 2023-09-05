# jointcal test data

This repository contains data to test the performance of the [jointcal](http://github.com/lsst/jointcal) product. Test data for jointcal needs to have multiple visits of the same field, processed with the LSST DM stack to get VisitInfo metadata and source catalogs.
See the note in the dependency file `ups/testdata_jointcal.table` regarding the obs package dependencies of this product: this package has no explicit dependencies, but does require the appropriate obs packages to be `setup` in order to use the contained butler repositories.

Individual sets of testing data should be placed in their own directories, as butler-accessible repositories. They can then be loaded in testing code via the butler and the catalogs fed to jointcal. The image HDUs in the calexps in each repo have been set to identically zero and then compressed to save space, since the calexps are only needed for their metadata.

The directories contained in this repository are listed below, with a description of their contents.

The `cfht`, and `hsc` directories each contain a `ref_cats/` directory with `gaia_dr2_20191105` and `ps1_pv3_3pi_20170110` reference catalogs (plus `sdss-dr9-fink-v5b` for `cfht`) in indexed HTM format.
The Gaia and PS1 refcats were extracted from the respective refcats on lsst-dev using the `scripts/extract-refcat-shards.py` script.
The sdss refcat was copied from the `validation_data_cfht/` repository's `ref_cats/` directory.

## cfht

Source catalogs, metadata, and zeroed+compressed images derived from [validation_data_cfht](https://github.com/lsst/validation_data_cfht), processed with weekly `w_2021_w46`.
These data contain two visits `(849375, 850587)`, centered at `(214.884832, 52.6622199)`, each containing detectors `(12,13,14,21,22,23)`.

The script `scripts/compress_jointcal_cfht_test_data.py` extracts the relevant files from the `validation_data_cfht` repository, removes the pixel-level data, and removes detectors that are not needed from the output and sourceTable_visit files.
`scripts/export_gen3_cfht.py` will produce the necessary gen3 `exports.yaml` file, to allow trivial importing of these output files into a butler repo for testing.

The parquet file columns were renamed according to the new schema (replacing leading capitals with leading lower case; replacing ``ccd`` with ``detector``; replacing ``filter`` with ``physical_filter`` and ``band``) in DM-31889, by running python scripts/rename_sourcetable_columns.py.

## hsc

Source catalogs, metadata, and zeroed+compressed images taken from the `w_2020_14`processing run of the HSC RC dataset, available at `lsst-dev:/datasets/hsc/repo/rerun/RC/w_2020_14/DM-24359-sfm`.
These data contain 10 visits, 5 in each of `HSC-R` and `HSC-I`, centered at `(337.710899, +0.807006)`, with these detectors for each visit (first 5 r, second 5 i):

```
    34648: [51, 59, 67],
    34690: [48, 56, 64],
    34714: [13, 18, 19],
    34674: [15, 21, 28],
    34670: [92, 93, 97],
    36140: [51, 59, 67],
    35892: [23, 31, 39],
    36192: [8, 14, 20],
    36260: [2, 7, 13],
    36236: [87, 93, 98]
```

The included `scripts/compress_jointcal_hsc_test_data.py` file copies the necessary data from the butler repo, removes pixel-level data from the calexps and compresses the source catalogs, and extracts the relevant portions of the sqlite3 registry.
The src catalogs have been converted to sourceTable_visit parquet tables for quick ingest.  The included `scripts/make_sourcetables_jointcal_hsc_test_data.sh` file runs the appropriate transformation tasks and cleans up temporary files.
The raw image files, which are necessary for gen3 conversion, are copied into the repo by `scripts/compress_jointcal_hsc_raw_test_data.py` which copies in the necessary data, and removes pixel-level data from the raws.

### gen3 conversion

To create a gen3 repo from the existing gen2 repo and make the exports file that is used to regenerate the gen3 repo for tests:

* Run `scripts/convert_gen2_to_gen3_hsc.sh` to create a gen3 repo from the existing gen2 repo (in the same location) and reduce the size of the created sky map.
* Run `scripts/export_gen3_hsc.py` to export the gen3 repo to `hsc/exports.yaml`.

We only commit `hsc/exports.yaml` and the "unbounded" `camera` and `transmission_*` calibrations to git, because those (and the in-place files in the gen2 repo) are all that is necessary to reconstruct a gen3 repo for testing.
Each test has to start with a fresh repo anyway, and creating one from an exports file is fast.

The parquet file columns were renamed according to the new schema (replacing leading capitals with leading lower case; replacing ``ccd`` with ``detector``; replacing ``filter`` with ``physical_filter`` and ``band``) in DM-31889, by running `python scripts/rename_sourcetable_columns.py`.

The skyMap was added to the export file (accidentally left behind) with DM-36726.

And the isolated star catalogs were added on DM-36726 by first rerunning `python scripts/rename_sourcetable_columns.py` and followed with `python scripts/make_hsc_star_associations.py` with some additional by-hand modification of the exports.yaml described in that script.
At the same time, the visitSummary tables were updated to include (blank) skyBg sky background values which are optionally used by fgcmcal, by running `python scripts/update_hsc_visit_summaries.py`.

## LATISS

On August 24, 2023, a subset of the dithered LATISS photometric survey data was extracted from the `LATISS/runs/AUXTEL_DRP_IMAGING_2023-08A-07AB-05AB/w_2023_33/PREOPS-3613` collection on `/repo/embargo`.
First, a new isolated star catalog needed to be created with only a subset of data, using `w_2023_33`:

* `pipetask run -j 1 -b /repo/embargo -i LATISS/runs/AUXTEL_DRP_IMAGING_2023-08A-07AB-05AB/w_2023_33/PREOPS-3613 -o u/erykoff/LATISS/2023-08A-07AB-05AB/testdata_subset -p $DRP_PIPE_DIR/pipelines/LATISS/DRP.yaml#isolatedStarAssociation -d "visit in (2023051100320, 2023051100357, 2023051100390, 2023051100406, 2023051100448, 2023051100454, 2023051100278, 2023051100473, 2023051100263, 2023051100509, 2023051100304, 2023051100431, 2023051100547, 2023051100379, 2023051100495, 2023051100489, 2023051100401, 2023051100280, 2023051100303, 2023051100508)"`

Then we run the following commands in the `testdata_jointcal` root directory:

* `mkdir -p latiss/testdata`
* `python scripts/export_latiss_from_embargo.py`

## Git LFS

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html)
explains how to set up Git LFS for LSST development.
