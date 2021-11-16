import lsst.pipe.tasks.postprocess
assert type(config)==lsst.pipe.tasks.postprocess.WriteSourceTableConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.postprocess.WriteSourceTableConfig' % (type(config).__module__, type(config).__name__)
# Add local photoCalib columns from the calexp.photoCalib? Should only set True if generating Source Tables from older src tables which do not already have local calib columns
config.doApplyExternalPhotoCalib=True

# Add local WCS columns from the calexp.wcs? Should only set True if generating Source Tables from older src tables which do not already have local calib columns
config.doApplyExternalSkyWcs=True

