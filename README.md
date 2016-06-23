jointcal validation data
========================

This repository contains data to validate the performance of the [jointcal](http://github.com/lsst/jointcal) product.

Individual sets of testing data should be placed in their own directories, as butler-accessible repositories. They can then be loaded in testing code via the butler and the catalogs fed to jointcal.

twinkles1
---------

Twinkles is a synthetic survey produced for the Dark Energy Science Collaboration's first data challenge. The intent with Twinkles is to test detection and analysis of type Ia SNe and strongly lensed quasars. Stars, galaxies and solar system objects were drawn from the LSST universe model. In addition to the standard universe, higher than natural rates of both type Ia SNe and lensed quasar systems were planted in the images. This repository includes one chip for the first 10 r-band visits in a deep drilling field (all observations are on the same night). See [here](https://github.com/DarkEnergyScienceCollaboration/Twinkles/blob/master/doc/Design.md) for a more complete description of the full survey.

twinkles1_and_index
-------------------

Astrometry index catalogs corresponding to the twinkles1 data.

Git LFS
-------

To clone and use this repository, you'll need Git Large File Storage (LFS).

Our [Developer Guide](http://developer.lsst.io/en/latest/tools/git_lfs.html)
explains how to set up Git LFS for LSST development.
