from lsst.obs.subaru import HyperSuprimeCam

testCollection = HyperSuprimeCam.makeCollectionName('testdata')
extraCollection = HyperSuprimeCam.makeCollectionName('extra')

config.refCats.append("gaia_dr2_20200414")
config.runs["gaia_dr2_20200414"] = "refcats"

config.runs["calexp"] = testCollection
config.runs["src"] = testCollection
config.runs["sourceTable_visit"] = testCollection
config.runs["visitSummary"] = testCollection
config.runs["src_schema"] = testCollection
config.runs["icSrc_schema"] = testCollection
config.runs["packages"] = extraCollection
config.runs["singleFrameDriver_config"] = extraCollection
config.runs["skyCorr_config"] = extraCollection
config.runs["transformSourceTable_config"] = extraCollection
config.runs["writeSourceTable_config"] = extraCollection

# Default HSC conversion config assumes there are brightObjectMask
# datasets landing in an HSC/masks collection, and then sets the
# option below to include HSC/masks in the HSC/defaults umbrella
# collection.  The repo we're converting here has no brightObjectMask
# datasets, so that collection won't exist.
config.extraUmbrellaChildren = []
