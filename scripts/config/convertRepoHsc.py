from lsst.obs.subaru import HyperSuprimeCam

testCollection = HyperSuprimeCam.makeCollectionName('testdata')
extraCollection = HyperSuprimeCam.makeCollectionName('extra')

config.runs["calexp"] = testCollection
config.runs["src"] = testCollection
config.runs["sourceTable_visit"] = testCollection
config.runs["src_schema"] = testCollection
config.runs["icSrc_schema"] = testCollection
config.runs["packages"] = extraCollection
config.runs["singleFrameDriver_config"] = extraCollection
config.runs["skyCorr_config"] = extraCollection
config.runs["transformSourceTable_config"] = extraCollection
config.runs["writeSourceTable_config"] = extraCollection
