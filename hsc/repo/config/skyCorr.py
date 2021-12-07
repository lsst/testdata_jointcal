import lsst.pipe.drivers.skyCorrection
assert type(config)==lsst.pipe.drivers.skyCorrection.SkyCorrectionConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.skyCorrection.SkyCorrectionConfig' % (type(config).__module__, type(config).__name__)
import lsst.pipe.base.config
import lsst.meas.algorithms.subtractBackground
import lsst.pipe.drivers.background
import lsst.meas.algorithms.detection
# Flag to enable/disable metadata saving for a task, enabled by default.
config.saveMetadata=True

# Bin size in x
config.bgModel.xSize=122.88

# Bin size in y
config.bgModel.ySize=122.88

# Pixel size in same units as xSize/ySize
config.bgModel.pixelSize=0.015

# Minimum fraction of bin size for good measurement
config.bgModel.minFrac=0.1

# Mask planes to treat as bad
config.bgModel.mask=['BAD', 'SAT', 'INTRP', 'DETECTED', 'DETECTED_NEGATIVE', 'EDGE', 'NO_DATA']

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.bgModel.interpolation='AKIMA_SPLINE'

# Do smoothing?
config.bgModel.doSmooth=False

# Smoothing scale, as a multiple of the bin size
config.bgModel.smoothScale=2.0

# Binning to use for CCD background model (pixels)
config.bgModel.binning=64

# Bin size in x
config.bgModel2.xSize=3.84

# Bin size in y
config.bgModel2.ySize=3.84

# Pixel size in same units as xSize/ySize
config.bgModel2.pixelSize=0.015

# Minimum fraction of bin size for good measurement
config.bgModel2.minFrac=0.5

# Mask planes to treat as bad
config.bgModel2.mask=['BAD', 'SAT', 'INTRP', 'DETECTED', 'DETECTED_NEGATIVE', 'EDGE', 'NO_DATA']

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.bgModel2.interpolation='AKIMA_SPLINE'

# Do smoothing?
config.bgModel2.doSmooth=True

# Smoothing scale, as a multiple of the bin size
config.bgModel2.smoothScale=1.0

# Binning to use for CCD background model (pixels)
config.bgModel2.binning=64

# k-sigma rejection iterations for sky scale
config.sky.skyIter=3

# k-sigma rejection threshold for sky scale
config.sky.skyRej=3.0

# type of statistic to use for grid points
config.sky.background.statistic='MEANCLIP'

# Superpixel size in x
config.sky.background.xBinSize=32

# Superpixel size in y
config.sky.background.yBinSize=32

# How to interpolate the background values. This maps to an enum; see afw::math::Background
config.sky.background.algorithm='NATURAL_SPLINE'

# Names of mask planes to ignore while estimating the background
config.sky.background.mask=['SAT', 'BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Number of samples in x for scaling sky frame
config.sky.xNumSamples=4

# Number of samples in y for scaling sky frame
config.sky.yNumSamples=4

# type of statistic to use for grid points
config.sky.stats.statistic='MEANCLIP'

# Clipping threshold for background
config.sky.stats.clip=3.0

# Clipping iterations for background
config.sky.stats.nIter=3

# Mask planes to reject
config.sky.stats.mask=['SAT', 'DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'NO_DATA']

# Number of iterations
config.maskObjects.nIter=3

# type of statistic to use for grid points
config.maskObjects.subtractBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.subtractBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.subtractBackground.binSize=1024

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.subtractBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.subtractBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.subtractBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.subtractBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.subtractBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.subtractBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.subtractBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.subtractBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.subtractBackground.weighting=True

# detected sources with fewer than the specified number of pixels will be ignored
config.maskObjects.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.maskObjects.detection.isotropicGrow=False

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.maskObjects.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.maskObjects.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.maskObjects.detection.returnOriginalFootprints=False

# Threshold for footprints; exact meaning and units depend on thresholdType.
config.maskObjects.detection.thresholdValue=2.5

# Include threshold relative to thresholdValue
config.maskObjects.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.maskObjects.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.maskObjects.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.maskObjects.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.maskObjects.detection.reEstimateBackground=False

# type of statistic to use for grid points
config.maskObjects.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.detection.background.weighting=True

# type of statistic to use for grid points
config.maskObjects.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.maskObjects.detection.doTempLocalBackground=False

# type of statistic to use for grid points
config.maskObjects.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.maskObjects.detection.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.maskObjects.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.maskObjects.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.maskObjects.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Detection threshold (standard deviations)
config.maskObjects.detectSigma=5.0

# Interpolate when removing objects?
config.maskObjects.doInterpolate=True

# type of statistic to use for grid points
config.maskObjects.interpolate.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.maskObjects.interpolate.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.maskObjects.interpolate.binSize=256

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.maskObjects.interpolate.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.maskObjects.interpolate.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.maskObjects.interpolate.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.maskObjects.interpolate.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.maskObjects.interpolate.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.maskObjects.interpolate.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.maskObjects.interpolate.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.maskObjects.interpolate.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.maskObjects.interpolate.weighting=True

# Mask objects to find good sky?
config.doMaskObjects=True

# Do background model subtraction?
config.doBgModel=True

# Do cleanup background model subtraction?
config.doBgModel2=True

# Do sky frame subtraction?
config.doSky=True

# Binning factor for constructing focal-plane images
config.binning=8

# Should be set to fakes_calexp if you want to process calexps with fakes in.
config.calexpType='calexp'

# name for connection rawLinker
config.connections.rawLinker='raw'

# name for connection calExpArray
config.connections.calExpArray='calexp'

# name for connection calBkgArray
config.connections.calBkgArray='calexpBackground'

# name for connection camera
config.connections.camera='camera'

# name for connection skyCalibs
config.connections.skyCalibs='sky'

# name for connection calExpCamera
config.connections.calExpCamera='calexp_camera'

# name for connection skyCorr
config.connections.skyCorr='skyCorr'

