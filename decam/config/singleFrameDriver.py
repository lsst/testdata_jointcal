import lsst.pipe.drivers.singleFrameDriver
assert type(config)==lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig, 'config is of type %s.%s instead of lsst.pipe.drivers.singleFrameDriver.SingleFrameDriverConfig' % (type(config).__module__, type(config).__name__)
import lsst.obs.decam.decamNullIsr
import lsst.meas.modelfit.cmodel.cmodelContinued
import lsst.meas.base.plugins
import lsst.meas.base.blendedness
import lsst.meas.algorithms.findCosmicRaysConfig
import lsst.pipe.tasks.photoCal
import lsst.shapelet.matrixBuilder
import lsst.meas.algorithms.installGaussianPsf
import lsst.meas.algorithms.skyObjects
import lsst.meas.modelfit.unitTransformedLikelihood
import lsst.pipe.tasks.postprocess
import lsst.meas.modelfit.optimizer
import lsst.meas.algorithms.detection
import lsst.meas.modelfit.cmodel.cmodel
import lsst.shapelet.radialProfile.radialProfile
import lsst.pipe.base.config
import lsst.shapelet.multiShapeletFunction.multiShapeletFunctionContinued
import lsst.meas.base.wrappers
import lsst.shapelet.version
import lsst.meas.algorithms.pcaPsfDeterminer
import lsst.pipe.tasks.characterizeImage
import lsst.meas.base.apertureFlux
import lsst.meas.modelfit.common
import lsst.shapelet.hermiteTransformMatrix
import lsst.shapelet
import lsst.meas.modelfit.adaptiveImportanceSampler
import lsst.meas.astrom.matchPessimisticB
import lsst.meas.astrom.directMatch
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegion
import lsst.meas.algorithms.gaussianPsfFactory
import lsst.shapelet.functorKeys
import lsst.meas.modelfit.model
import lsst.meas.base.classification
import lsst.shapelet.constants.constants
import lsst.meas.base.gaussianFlux
import lsst.meas.algorithms.measureApCorr
import lsst.meas.modelfit.version
import lsst.meas.base.sfm
import lsst.meas.base.localBackground
import lsst.meas.modelfit.psf.psfContinued
import lsst.afw.math.chebyshevBoundedFieldConfig
import lsst.pipe.tasks.fakes
import lsst.meas.modelfit.mixture
import lsst.meas.algorithms.reserveSourcesTask
import lsst.pipe.tasks.interpImage
import lsst.shapelet.basisEvaluator
import lsst.shapelet.radialProfile
import lsst.meas.modelfit.optimizer.optimizer
import lsst.meas.modelfit.psf
import lsst.meas.algorithms.matcherSourceSelector
import lsst.meas.algorithms.objectSizeStarSelector
import lsst.meas.astrom.astrometry
import lsst.meas.modelfit.truncatedGaussian
import lsst.meas.modelfit.pixelFitRegion
import lsst.meas.base.noiseReplacer
import lsst.shapelet.tractor
import lsst.shapelet.gaussHermiteConvolution
import lsst.meas.base.peakLikelihoodFlux
import lsst.meas.modelfit.unitSystem
import lsst.meas.modelfit.psf.psf
import lsst.meas.algorithms.loadIndexedReferenceObjects
import lsst.shapelet.multiShapeletBasis
import lsst.shapelet.multiShapeletFunction
import lsst.pipe.tasks.repair
import lsst.meas.modelfit
import lsst.meas.base.baseMeasurement
import lsst.shapelet.gaussHermiteProjection
import lsst.pipe.tasks.processCcd
import lsst.meas.modelfit.optimizer.optimizerContinued
import lsst.meas.base.sdssShape
import lsst.pipe.tasks.calibrate
import lsst.meas.base.scaledApertureFlux
import lsst.meas.modelfit.priors
import lsst.meas.modelfit.priors.priorsContinued
import lsst.shapelet.constants
import lsst.shapelet.constants.constantsContinued
import lsst.shapelet.shapeletFunction.shapeletFunctionContinued
import lsst.shapelet.shapeletFunction
import lsst.meas.algorithms.astrometrySourceSelector
import lsst.meas.modelfit.cmodel
import lsst.meas.base.pixelFlags
import lsst.meas.astrom.fitTanSipWcs
import lsst.obs.decam.crosstalk
import lsst.meas.base.applyApCorr
import lsst.meas.modelfit.likelihood
import lsst.meas.algorithms.makePsfCandidates
import lsst.pipe.tasks.colorterms
import lsst.meas.modelfit.multiModel
import lsst.meas.algorithms.sourceSelector
import lsst.shapelet.radialProfile.radialProfileContinued
import lsst.meas.algorithms.subtractBackground
import lsst.meas.astrom.ref_match
import lsst.meas.algorithms.flaggedSourceSelector
import lsst.meas.base.footprintArea
import lsst.meas.base.catalogCalculation
import lsst.meas.base.naiveCentroid
import lsst.meas.modelfit.pixelFitRegion.pixelFitRegionContinued
import lsst.meas.base.sdssCentroid
import lsst.shapelet.multiShapeletFunction.multiShapeletFunction
import lsst.meas.deblender.deblend
import lsst.pipe.tasks.measurePsf
import lsst.meas.modelfit.integrals
import lsst.meas.modelfit.sampler
import lsst.meas.modelfit.priors.priors
import lsst.shapelet.generator
import lsst.meas.base.psfFlux
import lsst.shapelet.shapeletFunction.shapeletFunction
config.processCcd.isr.retarget(target=lsst.obs.decam.decamNullIsr.DecamNullIsrTask, ConfigClass=lsst.obs.decam.decamNullIsr.DecamNullIsrConfig)

# Persist loaded data as a postISRCCD? The default is false, to avoid duplicating data.
config.processCcd.isr.doWrite=False

# Dataset type for input data; read by ProcessCcdTask; users will typically leave this alone
config.processCcd.isr.datasetType='instcal'

# Flag to enable/disable metadata saving for a task, enabled by default.
config.processCcd.charImage.saveMetadata=True

# Measure PSF? If False then for all subsequent operations use either existing PSF model when present, or install simple PSF model when not (see installSimplePsf config options)
config.processCcd.charImage.doMeasurePsf=True

# Persist results?
config.processCcd.charImage.doWrite=True

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.processCcd.charImage.doWriteExposure=False

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf is True then 2 should be plenty; otherwise more may be wanted.
config.processCcd.charImage.psfIterations=2

# type of statistic to use for grid points
config.processCcd.charImage.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.background.weighting=True

# detected sources with fewer than the specified number of pixels will be ignored
config.processCcd.charImage.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.charImage.detection.isotropicGrow=False

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.processCcd.charImage.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.processCcd.charImage.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.charImage.detection.returnOriginalFootprints=False

# Threshold for footprints; exact meaning and units depend on thresholdType.
config.processCcd.charImage.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.processCcd.charImage.detection.includeThresholdMultiplier=10.0

# specifies the desired flavor of Threshold
config.processCcd.charImage.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.processCcd.charImage.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.processCcd.charImage.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.processCcd.charImage.detection.reEstimateBackground=True

# type of statistic to use for grid points
config.processCcd.charImage.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.background.weighting=True

# type of statistic to use for grid points
config.processCcd.charImage.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.processCcd.charImage.detection.doTempLocalBackground=False

# type of statistic to use for grid points
config.processCcd.charImage.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.processCcd.charImage.detection.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.processCcd.charImage.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.processCcd.charImage.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.processCcd.charImage.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Run deblender input exposure
config.processCcd.charImage.doDeblend=False

# What to do when a peak to be deblended is close to the edge of the image
config.processCcd.charImage.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
config.processCcd.charImage.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.processCcd.charImage.deblend.assignStrayFlux=True

# How to split flux among peaks
config.processCcd.charImage.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.charImage.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.charImage.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.charImage.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.charImage.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.charImage.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.maxFootprintSize=0

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.charImage.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.processCcd.charImage.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.processCcd.charImage.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.processCcd.charImage.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.charImage.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.processCcd.charImage.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.charImage.deblend.maskLimits={}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.processCcd.charImage.deblend.weightTemplates=False

# Try to remove similar templates?
config.processCcd.charImage.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.processCcd.charImage.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.processCcd.charImage.deblend.medianSmoothTemplate=True

# the name of the centroiding algorithm used to set source x,y
config.processCcd.charImage.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.processCcd.charImage.measurement.slots.shape='base_SdssShape'

# the name of the algorithm used to set PSF moments parameters
config.processCcd.charImage.measurement.slots.psfShape='base_SdssShape_psf'

# the name of the algorithm used to set the source aperture instFlux slot
config.processCcd.charImage.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model instFlux slot
config.processCcd.charImage.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf instFlux slot
config.processCcd.charImage.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source Gaussian instFlux slot
config.processCcd.charImage.measurement.slots.gaussianFlux='base_GaussianFlux'

# the name of the instFlux measurement algorithm used for calibration
config.processCcd.charImage.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.processCcd.charImage.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.processCcd.charImage.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.processCcd.charImage.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation:
# >= 1: set the seed deterministically based on exposureId
# 0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.processCcd.charImage.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.processCcd.charImage.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.charImage.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.charImage.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.charImage.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.charImage.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.charImage.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.charImage.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.charImage.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.charImage.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.charImage.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.charImage.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.charImage.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.charImage.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.charImage.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_LocalPhotoCalib'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_LocalWcs'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.charImage.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.plugins['modelfit_CModel'].doMeasure=True

config.processCcd.charImage.measurement.plugins.names=['base_GaussianFlux', 'base_SdssShape', 'base_SdssCentroid', 'base_PsfFlux', 'base_CircularApertureFlux', 'base_PixelFlags']
# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.charImage.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.charImage.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.charImage.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.charImage.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.charImage.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.charImage.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.charImage.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.charImage.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.charImage.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.charImage.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.charImage.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.charImage.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_LocalPhotoCalib'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_LocalWcs'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.charImage.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.charImage.measurement.undeblended['modelfit_CModel'].doMeasure=True

config.processCcd.charImage.measurement.undeblended.names=[]
# Run subtasks to measure and apply aperture corrections
config.processCcd.charImage.doApCorr=True

# Field name prefix for the flux other measurements should be aperture corrected to match
config.processCcd.charImage.measureApCorr.refFluxName='slot_CalibFlux'

# Apply flux limit?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_instFlux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.measureApCorr.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.fluxField='base_PsfFlux_instFlux'

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['science'].signalToNoise.errField='base_PsfFlux_instFluxErr'

# Name of column for parent
config.processCcd.charImage.measureApCorr.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.measureApCorr.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.measureApCorr.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.measureApCorr.sourceSelector['references'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.measureApCorr.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.measureApCorr.sourceSelector['references'].colorLimits={}
# Apply flux limit to Psf Candidate selection?
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].doFluxLimit=True

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].fluxMax=0.0

# Apply signal-to-noise (i.e. flux/fluxErr) limit to Psf Candidate selection?
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].doSignalToNoiseLimit=False

# specify the minimum signal-to-noise for good Psf Candidates
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].signalToNoiseMin=20.0

# specify the maximum signal-to-noise for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].signalToNoiseMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_instFlux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.measureApCorr.sourceSelector['flagged'].field='calib_psf_used'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measureApCorr.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measureApCorr.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measureApCorr.sourceSelector['matcher'].minSnr=40.0

# Exclude objects that have saturated, interpolated, or edge pixels using PixelFlags. For matchOptimisticB set this to False to recover previous matcher selector behavior.
config.processCcd.charImage.measureApCorr.sourceSelector['matcher'].excludePixelFlags=True

config.processCcd.charImage.measureApCorr.sourceSelector.name='flagged'
# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
config.processCcd.charImage.measureApCorr.minDegreesOfFreedom=1

# maximum Chebyshev function order in x
config.processCcd.charImage.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.processCcd.charImage.measureApCorr.fitConfig.orderY=2

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.processCcd.charImage.measureApCorr.fitConfig.triangular=True

# Number of iterations for sigma clipping
config.processCcd.charImage.measureApCorr.numIter=4

# Number of standard devisations to clip at
config.processCcd.charImage.measureApCorr.numSigmaClip=3.0

# Allow these measurement algorithms to fail without an exception
config.processCcd.charImage.measureApCorr.allowFailure=[]

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.charImage.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.charImage.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.processCcd.charImage.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for modelFlux error
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.processCcd.charImage.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.processCcd.charImage.catalogCalculation.plugins.names=['base_FootprintArea', 'base_ClassificationExtendedness']
# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.processCcd.charImage.useSimplePsf=True

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.processCcd.charImage.installSimplePsf.fwhm=3.5322300675464238

# Width and height of PSF model, in pixels. Must be odd.
config.processCcd.charImage.installSimplePsf.width=11

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.charImage.refObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.charImage.refObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.charImage.refObjLoader.filterMap={'u': 'g', 'Y': 'y'}

# Require that the fields needed to correct proper motion (epoch, pm_ra and pm_dec) are present?
config.processCcd.charImage.refObjLoader.requireProperMotion=False

# Name of the ingested reference dataset
config.processCcd.charImage.refObjLoader.ref_dataset_name='ps1_pv3_3pi_20170110'

# Number of bright stars to use. Sets the max number of patterns that can be tested.
config.processCcd.charImage.ref_match.matcher.numBrightStars=200

# Minimum number of matched pairs; see also minFracMatchedPairs.
config.processCcd.charImage.ref_match.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs.
config.processCcd.charImage.ref_match.matcher.minFracMatchedPairs=0.3

# Number of softening iterations in matcher.
config.processCcd.charImage.ref_match.matcher.matcherIterations=5

# Maximum allowed shift of WCS, due to matching (pixel). When changing this value, the LoadReferenceObjectsConfig.pixelMargin should also be updated.
config.processCcd.charImage.ref_match.matcher.maxOffsetPix=300

# Rotation angle allowed between sources and position reference objects (degrees).
config.processCcd.charImage.ref_match.matcher.maxRotationDeg=1.0

# Number of points to define a shape for matching.
config.processCcd.charImage.ref_match.matcher.numPointsForShape=6

# Number of points to try for creating a shape. This value should be greater than or equal to numPointsForShape. Besides loosening the signal to noise cut in the 'matcher' SourceSelector, increasing this number will solve CCDs where no match was found.
config.processCcd.charImage.ref_match.matcher.numPointsForShapeAttempt=6

# Distance in units of pixels to always consider a source-reference pair a match. This prevents the astrometric fitter from over-fitting and removing stars that should be matched and allows for inclusion of new matches as the wcs improves.
config.processCcd.charImage.ref_match.matcher.minMatchDistPixels=1.0

# Number of implied shift/rotations from patterns that must agree before it a given shift/rotation is accepted. This is only used after the first softening iteration fails and if both the number of reference and source objects is greater than numBrightStars.
config.processCcd.charImage.ref_match.matcher.numPatternConsensus=3

# If the available reference objects exceeds this number, consensus/pessimistic mode will enforced regardless of the number of available sources. Below this optimistic mode (exit at first match rather than requiring numPatternConsensus to be matched) can be used. If more sources are required to match, decrease the signal to noise cut in the sourceSelector.
config.processCcd.charImage.ref_match.matcher.numRefRequireConsensus=1000

# Maximum number of reference objects to use for the matcher. The absolute maximum allowed for is 2 ** 16 for memory reasons.
config.processCcd.charImage.ref_match.matcher.maxRefObjects=65536

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
config.processCcd.charImage.ref_match.matchDistanceSigma=2.0

# Apply flux limit?
config.processCcd.charImage.ref_match.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.ref_match.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.ref_match.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_instFlux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['science'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.ref_match.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelector['science'].signalToNoise.fluxField='slot_CalibFlux_instFlux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.sourceSelector['science'].signalToNoise.errField='slot_CalibFlux_instFluxErr'

# Name of column for parent
config.processCcd.charImage.ref_match.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.ref_match.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.ref_match.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.sourceSelector['references'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.ref_match.sourceSelector['references'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.ref_match.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['references'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['references'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.ref_match.sourceSelector['references'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.ref_match.sourceSelector['references'].colorLimits={}
# Apply flux limit to Psf Candidate selection?
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].doFluxLimit=True

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].fluxMax=0.0

# Apply signal-to-noise (i.e. flux/fluxErr) limit to Psf Candidate selection?
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].doSignalToNoiseLimit=False

# specify the minimum signal-to-noise for good Psf Candidates
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].signalToNoiseMin=20.0

# specify the maximum signal-to-noise for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].signalToNoiseMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_instFlux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.ref_match.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.ref_match.sourceSelector['flagged'].field='calib_psf_used'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.ref_match.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.ref_match.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.ref_match.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.ref_match.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.ref_match.sourceSelector['matcher'].minSnr=40.0

# Exclude objects that have saturated, interpolated, or edge pixels using PixelFlags. For matchOptimisticB set this to False to recover previous matcher selector behavior.
config.processCcd.charImage.ref_match.sourceSelector['matcher'].excludePixelFlags=True

config.processCcd.charImage.ref_match.sourceSelector.name='science'
# Apply magnitude limit?
config.processCcd.charImage.ref_match.referenceSelector.doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.ref_match.referenceSelector.doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.ref_match.referenceSelector.doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.ref_match.referenceSelector.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.ref_match.referenceSelector.doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelector.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelector.magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.referenceSelector.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.ref_match.referenceSelector.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.ref_match.referenceSelector.flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelector.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelector.unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.ref_match.referenceSelector.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelector.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelector.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.ref_match.referenceSelector.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.referenceSelector.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.ref_match.referenceSelector.magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.ref_match.referenceSelector.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.ref_match.referenceSelector.magError.magErrField='mag_err'

config.processCcd.charImage.ref_match.referenceSelector.colorLimits={}
# Source flux type to use in source selection.
config.processCcd.charImage.ref_match.sourceFluxType='Calib'

# Apply flux limit?
config.processCcd.charImage.measurePsf.starSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measurePsf.starSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.charImage.measurePsf.starSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].fluxLimit.fluxField='slot_CalibFlux_instFlux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.measurePsf.starSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.fluxField='base_PsfFlux_instFlux'

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['science'].signalToNoise.errField='base_PsfFlux_instFluxErr'

# Name of column for parent
config.processCcd.charImage.measurePsf.starSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.charImage.measurePsf.starSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.charImage.measurePsf.starSelector['references'].doFlags=False

# Apply unresolved limitation?
config.processCcd.charImage.measurePsf.starSelector['references'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.charImage.measurePsf.starSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.charImage.measurePsf.starSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.charImage.measurePsf.starSelector['references'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.charImage.measurePsf.starSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.charImage.measurePsf.starSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.charImage.measurePsf.starSelector['references'].magError.magErrField='mag_err'

config.processCcd.charImage.measurePsf.starSelector['references'].colorLimits={}
# Apply flux limit to Psf Candidate selection?
config.processCcd.charImage.measurePsf.starSelector['objectSize'].doFluxLimit=True

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['objectSize'].fluxMax=0.0

# Apply signal-to-noise (i.e. flux/fluxErr) limit to Psf Candidate selection?
config.processCcd.charImage.measurePsf.starSelector['objectSize'].doSignalToNoiseLimit=False

# specify the minimum signal-to-noise for good Psf Candidates
config.processCcd.charImage.measurePsf.starSelector['objectSize'].signalToNoiseMin=20.0

# specify the maximum signal-to-noise for good Psf Candidates (ignored if == 0)
config.processCcd.charImage.measurePsf.starSelector['objectSize'].signalToNoiseMax=0.0

# minimum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.charImage.measurePsf.starSelector['objectSize'].sourceFluxField='base_GaussianFlux_instFlux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.charImage.measurePsf.starSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.charImage.measurePsf.starSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.charImage.measurePsf.starSelector['flagged'].field='calib_psf_used'

# List of flags which cause a source to be rejected as bad
config.processCcd.charImage.measurePsf.starSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measurePsf.starSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measurePsf.starSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.charImage.measurePsf.starSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.charImage.measurePsf.starSelector['matcher'].minSnr=40.0

# Exclude objects that have saturated, interpolated, or edge pixels using PixelFlags. For matchOptimisticB set this to False to recover previous matcher selector behavior.
config.processCcd.charImage.measurePsf.starSelector['matcher'].excludePixelFlags=True

config.processCcd.charImage.measurePsf.starSelector.name='objectSize'
# size of the kernel to create
config.processCcd.charImage.measurePsf.makePsfCandidates.kernelSize=21

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.makePsfCandidates.borderWidth=0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Minimum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# Maximum radius of the kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# Use non-linear fitter for spatial variation of Kernel
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# number of eigen components for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# specify spatial order for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# size of cell used to determine PSF (pixels, column direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# number of stars per psf cell for PSF kernel creation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of stars per psf Cell for spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Should each PSF candidate be given the same weight, independent of magnitude?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# number of iterations of PSF candidate star list
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=3

# tolerance of spatial fitting
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# floor for variance is lam*data
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=2.0

# Rejection threshold (stdev) for candidates based on spatial fit
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].spatialReject=3.0

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# Reject candidates that are blended?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# Mask blends in image?
config.processCcd.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

config.processCcd.charImage.measurePsf.psfDeterminer.name='pca'
# Fraction of candidates to reserve from fitting; none if <= 0
config.processCcd.charImage.measurePsf.reserve.fraction=0.0

# This number will be added to the exposure ID to set the random seed for reserving candidates
config.processCcd.charImage.measurePsf.reserve.seed=1

# Interpolate over defects? (ignored unless you provide a list of defects)
config.processCcd.charImage.repair.doInterpolate=True

# Find and mask out cosmic rays?
config.processCcd.charImage.repair.doCosmicRay=True

# maximum number of contaminated pixels
config.processCcd.charImage.repair.cosmicray.nCrPixelMax=100000

# CRs must be > this many sky-sig above sky
config.processCcd.charImage.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.processCcd.charImage.repair.cosmicray.min_DN=150.0

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.processCcd.charImage.repair.cosmicray.cond3_fac2=0.6

# number of times to look for contaminated pixels near known CR pixels
config.processCcd.charImage.repair.cosmicray.niteration=3

# Don't interpolate over CR pixels
config.processCcd.charImage.repair.cosmicray.keepCRs=False

# type of statistic to use for grid points
config.processCcd.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.charImage.repair.cosmicray.background.binSize=100000

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.charImage.repair.cosmicray.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.charImage.repair.cosmicray.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.charImage.repair.cosmicray.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.charImage.repair.cosmicray.background.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.charImage.repair.cosmicray.background.weighting=True

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.processCcd.charImage.repair.interp.modelPsf.size=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.sizeFactor=3.0

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.processCcd.charImage.repair.interp.modelPsf.maxSize=None

# Default FWHM of Gaussian model of core of star (pixels)
config.processCcd.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# Add a Gaussian to represent wings?
config.processCcd.charImage.repair.interp.modelPsf.addWing=True

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.processCcd.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Smoothly taper to the fallback value at the edge of the image?
config.processCcd.charImage.repair.interp.useFallbackValueAtEdge=True

# Type of statistic to calculate edge fallbackValue for interpolation
config.processCcd.charImage.repair.interp.fallbackValueType='MEANCLIP'

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.processCcd.charImage.repair.interp.fallbackUserValue=0.0

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.processCcd.charImage.repair.interp.negativeFallbackAllowed=True

# Transpose image before interpolating? This allows the interpolation to act over columns instead of rows.
config.processCcd.charImage.repair.interp.transpose=False

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.charImage.checkUnitsParseStrict='raise'

# name for connection exposure
config.processCcd.charImage.connections.exposure='postISRCCD'

# name for connection characterized
config.processCcd.charImage.connections.characterized='icExp'

# name for connection sourceCat
config.processCcd.charImage.connections.sourceCat='icSrc'

# name for connection backgroundModel
config.processCcd.charImage.connections.backgroundModel='icExpBackground'

# name for connection outputSchema
config.processCcd.charImage.connections.outputSchema='icSrc_schema'

# Perform calibration?
config.processCcd.doCalibrate=True

# Flag to enable/disable metadata saving for a task, enabled by default.
config.processCcd.calibrate.saveMetadata=True

# Save calibration results?
config.processCcd.calibrate.doWrite=True

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.processCcd.calibrate.doWriteHeavyFootprintsInSources=True

# Write reference matches (ignored if doWrite or doAstrometry false)?
config.processCcd.calibrate.doWriteMatches=True

# Write reference matches in denormalized format? This format uses more disk space, but is more convenient to read. Ignored if doWriteMatches=False or doWrite=False.
config.processCcd.calibrate.doWriteMatchesDenormalized=False

# Perform astrometric calibration?
config.processCcd.calibrate.doAstrometry=True

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.calibrate.astromRefObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.calibrate.astromRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.calibrate.astromRefObjLoader.filterMap={'u': 'phot_g_mean', 'g': 'phot_g_mean', 'r': 'phot_g_mean', 'i': 'phot_g_mean', 'z': 'phot_g_mean'}

# Require that the fields needed to correct proper motion (epoch, pm_ra and pm_dec) are present?
config.processCcd.calibrate.astromRefObjLoader.requireProperMotion=False

# Name of the ingested reference dataset
config.processCcd.calibrate.astromRefObjLoader.ref_dataset_name='gaia_dr2_20200414'

# Padding to add to 4 all edges of the bounding box (pixels)
config.processCcd.calibrate.photoRefObjLoader.pixelMargin=300

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.processCcd.calibrate.photoRefObjLoader.defaultFilter=''

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.processCcd.calibrate.photoRefObjLoader.filterMap={'u': 'g', 'Y': 'y'}

# Require that the fields needed to correct proper motion (epoch, pm_ra and pm_dec) are present?
config.processCcd.calibrate.photoRefObjLoader.requireProperMotion=False

# Name of the ingested reference dataset
config.processCcd.calibrate.photoRefObjLoader.ref_dataset_name='ps1_pv3_3pi_20170110'

# Number of bright stars to use. Sets the max number of patterns that can be tested.
config.processCcd.calibrate.astrometry.matcher.numBrightStars=200

# Minimum number of matched pairs; see also minFracMatchedPairs.
config.processCcd.calibrate.astrometry.matcher.minMatchedPairs=30

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs.
config.processCcd.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# Number of softening iterations in matcher.
config.processCcd.calibrate.astrometry.matcher.matcherIterations=5

# Maximum allowed shift of WCS, due to matching (pixel). When changing this value, the LoadReferenceObjectsConfig.pixelMargin should also be updated.
config.processCcd.calibrate.astrometry.matcher.maxOffsetPix=300

# Rotation angle allowed between sources and position reference objects (degrees).
config.processCcd.calibrate.astrometry.matcher.maxRotationDeg=1.0

# Number of points to define a shape for matching.
config.processCcd.calibrate.astrometry.matcher.numPointsForShape=6

# Number of points to try for creating a shape. This value should be greater than or equal to numPointsForShape. Besides loosening the signal to noise cut in the 'matcher' SourceSelector, increasing this number will solve CCDs where no match was found.
config.processCcd.calibrate.astrometry.matcher.numPointsForShapeAttempt=6

# Distance in units of pixels to always consider a source-reference pair a match. This prevents the astrometric fitter from over-fitting and removing stars that should be matched and allows for inclusion of new matches as the wcs improves.
config.processCcd.calibrate.astrometry.matcher.minMatchDistPixels=1.0

# Number of implied shift/rotations from patterns that must agree before it a given shift/rotation is accepted. This is only used after the first softening iteration fails and if both the number of reference and source objects is greater than numBrightStars.
config.processCcd.calibrate.astrometry.matcher.numPatternConsensus=3

# If the available reference objects exceeds this number, consensus/pessimistic mode will enforced regardless of the number of available sources. Below this optimistic mode (exit at first match rather than requiring numPatternConsensus to be matched) can be used. If more sources are required to match, decrease the signal to noise cut in the sourceSelector.
config.processCcd.calibrate.astrometry.matcher.numRefRequireConsensus=1000

# Maximum number of reference objects to use for the matcher. The absolute maximum allowed for is 2 ** 16 for memory reasons.
config.processCcd.calibrate.astrometry.matcher.maxRefObjects=65536

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
config.processCcd.calibrate.astrometry.matchDistanceSigma=2.0

# Apply flux limit?
config.processCcd.calibrate.astrometry.sourceSelector['science'].doFluxLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.sourceSelector['science'].doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.astrometry.sourceSelector['science'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.sourceSelector['science'].doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.calibrate.astrometry.sourceSelector['science'].doIsolated=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelector['science'].fluxLimit.fluxField='slot_CalibFlux_instFlux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelector['science'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelector['science'].flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_saturated', 'base_PsfFlux_flags']

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.calibrate.astrometry.sourceSelector['science'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['science'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelector['science'].signalToNoise.fluxField='base_PsfFlux_instFlux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.sourceSelector['science'].signalToNoise.errField='base_PsfFlux_instFluxErr'

# Name of column for parent
config.processCcd.calibrate.astrometry.sourceSelector['science'].isolated.parentName='parent'

# Name of column for nChild
config.processCcd.calibrate.astrometry.sourceSelector['science'].isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.calibrate.astrometry.sourceSelector['references'].doMagLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.sourceSelector['references'].doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.astrometry.sourceSelector['references'].doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.sourceSelector['references'].doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.astrometry.sourceSelector['references'].doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelector['references'].magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelector['references'].flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.sourceSelector['references'].flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.calibrate.astrometry.sourceSelector['references'].unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.sourceSelector['references'].signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.sourceSelector['references'].signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.sourceSelector['references'].magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.sourceSelector['references'].magError.magErrField='mag_err'

config.processCcd.calibrate.astrometry.sourceSelector['references'].colorLimits={}
# Apply flux limit to Psf Candidate selection?
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].doFluxLimit=True

# specify the minimum psfFlux for good Psf Candidates
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].fluxMin=12500.0

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].fluxMax=0.0

# Apply signal-to-noise (i.e. flux/fluxErr) limit to Psf Candidate selection?
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].doSignalToNoiseLimit=False

# specify the minimum signal-to-noise for good Psf Candidates
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].signalToNoiseMin=20.0

# specify the maximum signal-to-noise for good Psf Candidates (ignored if == 0)
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].signalToNoiseMax=0.0

# minimum width to include in histogram
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].widthMin=0.0

# maximum width to include in histogram
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].widthMax=10.0

# Name of field in Source to use for flux measurement
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].sourceFluxField='base_GaussianFlux_instFlux'

# Standard deviation of width allowed to be interpreted as good stars
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].widthStdAllowed=0.15

# Keep objects within this many sigma of cluster 0's median
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].nSigmaClip=2.0

# List of flags which cause a source to be rejected as bad
config.processCcd.calibrate.astrometry.sourceSelector['objectSize'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# Name of a flag field that is True for Sources that should be used.
config.processCcd.calibrate.astrometry.sourceSelector['flagged'].field='calib_psf_used'

# List of flags which cause a source to be rejected as bad
config.processCcd.calibrate.astrometry.sourceSelector['astrometry'].badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad']

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.sourceSelector['astrometry'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.calibrate.astrometry.sourceSelector['astrometry'].minSnr=10.0

# Type of source flux; typically one of Ap or Psf
config.processCcd.calibrate.astrometry.sourceSelector['matcher'].sourceFluxType='Ap'

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <= 0 for no limit
config.processCcd.calibrate.astrometry.sourceSelector['matcher'].minSnr=40.0

# Exclude objects that have saturated, interpolated, or edge pixels using PixelFlags. For matchOptimisticB set this to False to recover previous matcher selector behavior.
config.processCcd.calibrate.astrometry.sourceSelector['matcher'].excludePixelFlags=True

config.processCcd.calibrate.astrometry.sourceSelector.name='matcher'
# Apply magnitude limit?
config.processCcd.calibrate.astrometry.referenceSelector.doMagLimit=False

# Apply flag limitation?
config.processCcd.calibrate.astrometry.referenceSelector.doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.astrometry.referenceSelector.doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.astrometry.referenceSelector.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.astrometry.referenceSelector.doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelector.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelector.magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.referenceSelector.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.astrometry.referenceSelector.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.astrometry.referenceSelector.flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelector.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelector.unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.calibrate.astrometry.referenceSelector.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelector.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelector.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.astrometry.referenceSelector.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.referenceSelector.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.astrometry.referenceSelector.magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.astrometry.referenceSelector.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.astrometry.referenceSelector.magError.magErrField='mag_err'

config.processCcd.calibrate.astrometry.referenceSelector.colorLimits={}
# Source flux type to use in source selection.
config.processCcd.calibrate.astrometry.sourceFluxType='Ap'

# order of SIP polynomial
config.processCcd.calibrate.astrometry.wcsFitter.order=4

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
config.processCcd.calibrate.astrometry.wcsFitter.numIter=3

# number of rejection iterations
config.processCcd.calibrate.astrometry.wcsFitter.numRejIter=1

# Number of standard deviations for clipping level
config.processCcd.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
config.processCcd.calibrate.astrometry.wcsFitter.maxScatterArcsec=20.0

# If True then load reference objects and match sources but do not fit a WCS; this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.processCcd.calibrate.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
config.processCcd.calibrate.astrometry.maxIter=3

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
config.processCcd.calibrate.astrometry.minMatchDistanceArcSec=0.001

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.processCcd.calibrate.requireAstrometry=True

# Perform phometric calibration?
config.processCcd.calibrate.doPhotoCal=True

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.processCcd.calibrate.requirePhotoCal=True

# Matching radius, arcsec
config.processCcd.calibrate.photoCal.match.matchRadius=0.25

# Apply flux limit?
config.processCcd.calibrate.photoCal.match.sourceSelection.doFluxLimit=False

# Apply flag limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doFlags=True

# Apply unresolved limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doUnresolved=True

# Apply signal-to-noise limit?
config.processCcd.calibrate.photoCal.match.sourceSelection.doSignalToNoise=False

# Apply isolated limitation?
config.processCcd.calibrate.photoCal.match.sourceSelection.doIsolated=False

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.fluxLimit.fluxField='slot_CalibFlux_instFlux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.photoCal.match.sourceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.photoCal.match.sourceSelection.flags.bad=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.calibrate.photoCal.match.sourceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.fluxField='base_PsfFlux_instFlux'

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.sourceSelection.signalToNoise.errField='base_PsfFlux_instFluxErr'

# Name of column for parent
config.processCcd.calibrate.photoCal.match.sourceSelection.isolated.parentName='parent'

# Name of column for nChild
config.processCcd.calibrate.photoCal.match.sourceSelection.isolated.nChildName='deblend_nChild'

# Apply magnitude limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doMagLimit=False

# Apply flag limitation?
config.processCcd.calibrate.photoCal.match.referenceSelection.doFlags=False

# Apply unresolved limitation?
config.processCcd.calibrate.photoCal.match.referenceSelection.doUnresolved=False

# Apply signal-to-noise limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doSignalToNoise=False

# Apply magnitude error limit?
config.processCcd.calibrate.photoCal.match.referenceSelection.doMagError=False

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.magLimit.fluxField='flux'

# List of source flag fields that must be set for a source to be used.
config.processCcd.calibrate.photoCal.match.referenceSelection.flags.good=[]

# List of source flag fields that must NOT be set for a source to be used.
config.processCcd.calibrate.photoCal.match.referenceSelection.flags.bad=[]

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.unresolved.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.unresolved.maximum=0.5

# Name of column for star/galaxy separation
config.processCcd.calibrate.photoCal.match.referenceSelection.unresolved.name='base_ClassificationExtendedness_value'

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.maximum=None

# Name of the source flux field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.fluxField='flux'

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.signalToNoise.errField='flux_err'

# Select objects with value greater than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.minimum=None

# Select objects with value less than this
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.maximum=None

# Name of the source flux error field to use.
config.processCcd.calibrate.photoCal.match.referenceSelection.magError.magErrField='mag_err'

config.processCcd.calibrate.photoCal.match.referenceSelection.colorLimits={}
# Fraction of candidates to reserve from fitting; none if <= 0
config.processCcd.calibrate.photoCal.reserve.fraction=0.0

# This number will be added to the exposure ID to set the random seed for reserving candidates
config.processCcd.calibrate.photoCal.reserve.seed=1

# Name of the source instFlux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.processCcd.calibrate.photoCal.fluxField='slot_CalibFlux_instFlux'

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.processCcd.calibrate.photoCal.applyColorTerms=None

# maximum sigma to use when clipping
config.processCcd.calibrate.photoCal.sigmaMax=0.25

# clip at nSigma
config.processCcd.calibrate.photoCal.nSigma=3.0

# use median instead of mean to compute zeropoint
config.processCcd.calibrate.photoCal.useMedian=True

# number of iterations
config.processCcd.calibrate.photoCal.nIter=20

config.processCcd.calibrate.photoCal.colorterms.data={}
# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.processCcd.calibrate.photoCal.photoCatName='ps1_pv3_3pi_20170110'

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
config.processCcd.calibrate.photoCal.magErrFloor=0.0

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. Ignored if icSourceCat is not provided.
config.processCcd.calibrate.icSourceFieldsToCopy=['calib_psf_candidate', 'calib_psf_used', 'calib_psf_reserved']

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.processCcd.calibrate.matchRadiusPix=3.0

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.processCcd.calibrate.checkUnitsParseStrict='raise'

# detected sources with fewer than the specified number of pixels will be ignored
config.processCcd.calibrate.detection.minPixels=1

# Pixels should be grown as isotropically as possible (slower)
config.processCcd.calibrate.detection.isotropicGrow=False

# Grow all footprints at the same time? This allows disconnected footprints to merge.
config.processCcd.calibrate.detection.combinedGrow=True

# Grow detections by nSigmaToGrow * [PSF RMS width]; if 0 then do not grow
config.processCcd.calibrate.detection.nSigmaToGrow=2.4

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.processCcd.calibrate.detection.returnOriginalFootprints=False

# Threshold for footprints; exact meaning and units depend on thresholdType.
config.processCcd.calibrate.detection.thresholdValue=5.0

# Include threshold relative to thresholdValue
config.processCcd.calibrate.detection.includeThresholdMultiplier=1.0

# specifies the desired flavor of Threshold
config.processCcd.calibrate.detection.thresholdType='stdev'

# specifies whether to detect positive, or negative sources, or both
config.processCcd.calibrate.detection.thresholdPolarity='positive'

# Fiddle factor to add to the background; debugging only
config.processCcd.calibrate.detection.adjustBackground=0.0

# Estimate the background again after final source detection?
config.processCcd.calibrate.detection.reEstimateBackground=True

# type of statistic to use for grid points
config.processCcd.calibrate.detection.background.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.background.binSize=128

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.background.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.background.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.background.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.background.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.background.useApprox=True

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.background.weighting=True

# type of statistic to use for grid points
config.processCcd.calibrate.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.tempLocalBackground.binSize=64

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempLocalBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempLocalBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.tempLocalBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.tempLocalBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempLocalBackground.weighting=True

# Enable temporary local background subtraction? (see tempLocalBackground)
config.processCcd.calibrate.detection.doTempLocalBackground=False

# type of statistic to use for grid points
config.processCcd.calibrate.detection.tempWideBackground.statisticsProperty='MEANCLIP'

# behaviour if there are too few points in grid for requested interpolation style
config.processCcd.calibrate.detection.tempWideBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how large a region of the sky should be used for each background point
config.processCcd.calibrate.detection.tempWideBackground.binSize=512

# Sky region size to be used for each background point in X direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempWideBackground.binSizeX=0

# Sky region size to be used for each background point in Y direction. If 0, the binSize config is used.
config.processCcd.calibrate.detection.tempWideBackground.binSizeY=0

# how to interpolate the background values. This maps to an enum; see afw::math::Background
config.processCcd.calibrate.detection.tempWideBackground.algorithm='AKIMA_SPLINE'

# Names of mask planes to ignore while estimating the background
config.processCcd.calibrate.detection.tempWideBackground.ignoredPixelMask=['BAD', 'EDGE', 'NO_DATA']

# Ignore NaNs when estimating the background
config.processCcd.calibrate.detection.tempWideBackground.isNanSafe=False

# Use Approximate (Chebyshev) to model background.
config.processCcd.calibrate.detection.tempWideBackground.useApprox=False

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.approxOrderY=-1

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.processCcd.calibrate.detection.tempWideBackground.weighting=True

# Do temporary wide (large-scale) background subtraction before footprint detection?
config.processCcd.calibrate.detection.doTempWideBackground=False

# The maximum number of peaks in a Footprint before trying to replace its peaks using the temporary local background
config.processCcd.calibrate.detection.nPeaksMaxSimple=1

# Multiple of PSF RMS size to use for convolution kernel bounding box size; note that this is not a half-size. The size will be rounded up to the nearest odd integer
config.processCcd.calibrate.detection.nSigmaForKernel=7.0

# Mask planes to ignore when calculating statistics of image (for thresholdType=stdev)
config.processCcd.calibrate.detection.statsMask=['BAD', 'SAT', 'EDGE', 'NO_DATA']

# Run deblender input exposure
config.processCcd.calibrate.doDeblend=True

# What to do when a peak to be deblended is close to the edge of the image
config.processCcd.calibrate.deblend.edgeHandling='ramp'

# When the deblender should attribute stray flux to point sources
config.processCcd.calibrate.deblend.strayFluxToPointSources='necessary'

# Assign stray flux (not claimed by any child in the deblender) to deblend children.
config.processCcd.calibrate.deblend.assignStrayFlux=True

# How to split flux among peaks
config.processCcd.calibrate.deblend.strayFluxRule='trim'

# When splitting stray flux, clip fractions below this value to zero.
config.processCcd.calibrate.deblend.clipStrayFluxFraction=0.001

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq1=1.5

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.processCcd.calibrate.deblend.psfChisq2=1.5

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.processCcd.calibrate.deblend.psfChisq2b=1.5

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.processCcd.calibrate.deblend.maxNumberOfPeaks=0

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintArea=1000000

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.maxFootprintSize=2000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.processCcd.calibrate.deblend.minFootprintAxisRatio=0.0

# Mask name for footprints not deblended, or None
config.processCcd.calibrate.deblend.notDeblendedMask='NOT_DEBLENDED'

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
config.processCcd.calibrate.deblend.tinyFootprintSize=2

# Guarantee that all peaks produce a child source.
config.processCcd.calibrate.deblend.propagateAllPeaks=False

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.processCcd.calibrate.deblend.catchFailures=False

# Mask planes to ignore when performing statistics
config.processCcd.calibrate.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.processCcd.calibrate.deblend.maskLimits={}

# If true, a least-squares fit of the templates will be done to the full image. The templates will be re-weighted based on this fit.
config.processCcd.calibrate.deblend.weightTemplates=False

# Try to remove similar templates?
config.processCcd.calibrate.deblend.removeDegenerateTemplates=False

# If the dot product between two templates is larger than this value, we consider them to be describing the same object (i.e. they are degenerate).  If one of the objects has been labeled as a PSF it will be removed, otherwise the template with the lowest value will be removed.
config.processCcd.calibrate.deblend.maxTempDotProd=0.5

# Apply a smoothing filter to all of the template images
config.processCcd.calibrate.deblend.medianSmoothTemplate=True

# Generate sky sources?
config.processCcd.calibrate.doSkySources=True

# Avoid pixels masked with these mask planes
config.processCcd.calibrate.skySources.avoidMask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'NO_DATA']

# Number of pixels to grow the masked pixels when adding sky objects
config.processCcd.calibrate.skySources.growMask=0

# Radius, in pixels, of sky objects
config.processCcd.calibrate.skySources.sourceRadius=8.0

# Try to add this many sky objects
config.processCcd.calibrate.skySources.nSources=100

# Maximum number of trial sky object positions
# (default: nSkySources*nTrialSkySourcesMultiplier)
config.processCcd.calibrate.skySources.nTrialSources=None

# Set nTrialSkySources to
#     nSkySources*nTrialSkySourcesMultiplier
# if nTrialSkySources is None
config.processCcd.calibrate.skySources.nTrialSourcesMultiplier=5

# the name of the centroiding algorithm used to set source x,y
config.processCcd.calibrate.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set source moments parameters
config.processCcd.calibrate.measurement.slots.shape='base_SdssShape'

# the name of the algorithm used to set PSF moments parameters
config.processCcd.calibrate.measurement.slots.psfShape='base_SdssShape_psf'

# the name of the algorithm used to set the source aperture instFlux slot
config.processCcd.calibrate.measurement.slots.apFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source model instFlux slot
config.processCcd.calibrate.measurement.slots.modelFlux='modelfit_CModel'

# the name of the algorithm used to set the source psf instFlux slot
config.processCcd.calibrate.measurement.slots.psfFlux='base_PsfFlux'

# the name of the algorithm used to set the source Gaussian instFlux slot
config.processCcd.calibrate.measurement.slots.gaussianFlux='base_GaussianFlux'

# the name of the instFlux measurement algorithm used for calibration
config.processCcd.calibrate.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# When measuring, replace other detected footprints with noise?
config.processCcd.calibrate.measurement.doReplaceWithNoise=True

# How to choose mean and variance of the Gaussian noise we generate?
config.processCcd.calibrate.measurement.noiseReplacer.noiseSource='measure'

# Add ann offset to the generated noise.
config.processCcd.calibrate.measurement.noiseReplacer.noiseOffset=0.0

# The seed multiplier value to use for random number generation:
# >= 1: set the seed deterministically based on exposureId
# 0: fall back to the afw.math.Random default constructor (which uses a seed value of 1)
config.processCcd.calibrate.measurement.noiseReplacer.noiseSeedMultiplier=1

# Prefix to give undeblended plugins
config.processCcd.calibrate.measurement.undeblendedPrefix='undeblended_'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.calibrate.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.calibrate.measurement.plugins['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.plugins['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.calibrate.measurement.plugins['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.calibrate.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.calibrate.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.calibrate.measurement.plugins['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.calibrate.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.calibrate.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.calibrate.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_LocalPhotoCalib'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_LocalWcs'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.calibrate.measurement.plugins['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.plugins['modelfit_CModel'].doMeasure=True

config.processCcd.calibrate.measurement.plugins.names=['base_PixelFlags', 'modelfit_DoubleShapeletPsfApprox', 'base_GaussianFlux', 'base_Variance', 'base_LocalBackground', 'base_SdssShape', 'base_SdssCentroid', 'base_LocalWcs', 'base_CircularApertureFlux', 'base_PsfFlux', 'base_Blendedness', 'base_SkyCoord', 'base_LocalPhotoCalib', 'modelfit_CModel', 'base_NaiveCentroid']
# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.processCcd.calibrate.measurement.undeblended['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.processCcd.calibrate.measurement.undeblended['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.processCcd.calibrate.measurement.undeblended['base_GaussianFlux'].background=0.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].background=0.0

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.undeblended['base_NaiveCentroid'].maxDistToPeak=-1.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].doMeasure=True

# maximum allowed binning
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].binmax=16

# Do check that the centroid is contained in footprint.
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].doFootprintCheck=True

# If set > 0, Centroid Check also checks distance from footprint peak.
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].maxDistToPeak=-1.0

# if the peak's less than this insist on binning at least once
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.processCcd.calibrate.measurement.undeblended['base_SdssCentroid'].wfac=1.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].masksFpAnywhere=[]

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.processCcd.calibrate.measurement.undeblended['base_PixelFlags'].masksFpCenter=[]

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].background=0.0

# Whether to also compute the shape of the PSF model
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].doMeasurePsf=True

# Maximum number of iterations
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].maxIter=100

# Maximum centroid shift, limited to 2-10
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].maxShift=0.0

# Convergence tolerance for e1,e2
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].tol1=9.999999747378752e-06

# Convergence tolerance for FWHM
config.processCcd.calibrate.measurement.undeblended['base_SdssShape'].tol2=9.999999747378752e-05

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].doMeasure=True

# Scaling factor of PSF FWHM for aperture radius.
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].scale=3.14

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.undeblended['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].doMeasure=True

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].maxSincRadius=10.0

# Radius (in pixels) of apertures.
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.processCcd.calibrate.measurement.undeblended['base_CircularApertureFlux'].shiftKernel='lanczos5'

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doFlux=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].doShape=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.processCcd.calibrate.measurement.undeblended['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].doMeasure=True

# Inner radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].annulusInner=7.0

# Outer radius for background annulus as a multiple of the PSF sigma
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].annulusOuter=15.0

# Mask planes that indicate pixels that should be excluded from the measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].badMaskPlanes=['BAD', 'SAT', 'NO_DATA']

# Number of sigma-clipping iterations for background measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].bgIter=3

# Rejection threshold (in standard deviations) for background measurement
config.processCcd.calibrate.measurement.undeblended['base_LocalBackground'].bgRej=3.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.processCcd.calibrate.measurement.undeblended['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.processCcd.calibrate.measurement.undeblended['base_Variance'].scale=5.0

# Mask planes to ignore
config.processCcd.calibrate.measurement.undeblended['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_LocalPhotoCalib'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_LocalWcs'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].doMeasure=True

# Shapelet order of inner expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].innerOrder=2

# Don't allow the semi-major radius of any component to go above this fraction of the PSF image width
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].maxRadiusBoxFraction=0.4

# Don't allow the semi-minor radius of any component to drop below this value (pixels)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadius=1.0

# Don't allow the determinant radii of the two components to differ by less than this (pixels)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].minRadiusDiff=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].optimizer.trustRegionSolverTolerance=1e-08

# Shapelet order of outer expansion (0 == Gaussian)
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].outerOrder=1

# Initial outer Gaussian peak height divided by inner Gaussian peak height
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].peakRatio=0.1

# Initial outer radius divided by inner radius
config.processCcd.calibrate.measurement.undeblended['modelfit_DoubleShapeletPsfApprox'].radiusRatio=2.0

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models={}
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['SingleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleGaussian'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.order=-1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.order=2

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.order=1

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['DoubleShapelet'].wings.radiusPriorSigma=0.5

config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full']=lsst.meas.modelfit.GeneralPsfFitterConfig()
# Default value for the noiseSigma parameter in GeneralPsfFitter.apply()
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].defaultNoiseSigma=0.001

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusFactor=0.5

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].inner.radiusPriorSigma=0.5

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].optimizer.trustRegionSolverTolerance=1e-08

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.order=0

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusFactor=4.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].outer.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusFactor=1.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].primary.radiusPriorSigma=0.5

# sigma in an isotropic 2-d Gaussian prior on the conformal-shear ellipticity eta
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.ellipticityPriorSigma=0.3

# shapelet order for this component; negative to disable this component completely
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.order=4

# sigma (in pixels) in an isotropic 2-d Gaussian prior on the center of this shapelet component, relative to the center of the PSF image
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.positionPriorSigma=0.1

# Sets the fiducial radius of this component relative to the 'primary radius' of the PSF: either the second-moments radius of the PSF image (in an initial fit), or the radius of the primary component in a previous fit.  Ignored if the previous fit included this component (as then we can just use that radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusFactor=2.0

# sigma in a Gaussian prior on ln(radius/fiducialRadius)
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].models['Full'].wings.radiusPriorSigma=0.5

# a sequence of model names indicating which models should be fit, and their order
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].sequence=['DoubleShapelet']

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_GeneralShapeletPsfApprox'].doMeasure=True

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.nComponents=8

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.profileName='luv'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].dev.weightsMultiplier=1.0

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.nComponents=6

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.gradientThreshold=1e-05

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.maxOuterIterations=250

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.minTrustRadiusThreshold=1e-05

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.usePixelWeights=False

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].exp.weightsMultiplier=1.0

# If the 2nd-moments shape used to initialize the fit failed, use the PSF moments multiplied by this.  If <= 0.0, abort the fit early instead.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].fallbackInitialMomentsPsfFactor=1.5

# Whether to record the steps the optimizer takes (or just the number, if running as a plugin)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.doRecordHistory=True

# Whether to record the time spent in this stage
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.doRecordTime=True

# Softened core width for ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticityCore=0.001

# Width of exponential ellipticity distribution (conformal shear units).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.ellipticitySigma=0.3

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinInner=-6.0

# Minimum ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMinOuter=-6.001

# Mean of the Student's T distribution used for ln(radius) at large radius, and the transition point between a flat distribution and the Student's T.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusMu=-1.0

# Number of degrees of freedom for the Student's T distribution on ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusNu=50.0

# Width of the Student's T distribution in ln(radius).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.empiricalPriorConfig.logRadiusSigma=0.45

# Ellipticity magnitude (conformal shear units) at which the softened cutoff begins
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxInner=2.0

# Maximum ellipticity magnitude (conformal shear units)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.ellipticityMaxOuter=2.001

# ln(radius) at which the softened cutoff begins towards the maximum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxInner=3.0

# Maximum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMaxOuter=3.001

# ln(radius) at which the softened cutoff begins towards the minimum
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinInner=-6.0

# The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinMaxRatio=1.0

# Minimum ln(radius)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.linearPriorConfig.logRadiusMinOuter=-6.001

# Maximum radius used in approximating profile with Gaussians (0=default for this profile)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.maxRadius=0

# Number of Gaussian used to approximate the profile
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.nComponents=3

# whether to save all iterations for debugging purposes
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.doSaveIterations=False

# If the maximum of the gradient falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.gradientThreshold=0.001

# maximum number of iterations (i.e. function evaluations and trust region subproblems) per step
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxInnerIterations=20

# maximum number of steps
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.maxOuterIterations=500

# If the trust radius falls below this threshold, consider the algorithm converged
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.minTrustRadiusThreshold=0.01

# If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.noSR1Term=False

# absolute step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffAbsStep=0.0

# relative step size used for numerical derivatives (added to other steps)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffRelStep=0.0

# step size (in units of trust radius) used for numerical derivatives (added to relative step)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.numDiffTrustRadiusStep=0.1

# Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.skipSR1UpdateThreshold=1e-08

# steps with reduction ratio greater than this are accepted
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.stepAcceptThreshold=0.0

# when increase the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowFactor=2.0

# steps with reduction radio greater than this may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowReductionRatio=0.75

# steps with length this fraction of the trust radius may increase the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionGrowStepFraction=0.8

# the initial trust region will be set to this value
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionInitialSize=1.0

# when reducing the trust region size, multiply the radius by this factor
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkFactor=0.3333333333333333

# steps with reduction radio less than this will decrease the trust radius
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionShrinkReductionRatio=0.25

# value passed as the tolerance to solveTrustRegion
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.optimizer.trustRegionSolverTolerance=1e-08

# Name of the Prior that defines the model to fit (a filename in $MEAS_MODELFIT_DIR/data, with no extension), if priorSource='FILE'.  Ignored for forced fitting.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.priorName=''

# One of 'FILE', 'LINEAR', 'EMPIRICAL', or 'NONE', indicating whether the prior should be loaded from disk, created from one of the nested prior config/control objects, or None
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.priorSource='EMPIRICAL'

# Name of the shapelet.RadialProfile that defines the model to fit
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.profileName='lux'

# Use per-pixel variances as weights in the nonlinear fit (the final linear fit for flux never uses per-pixel variances)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.usePixelWeights=True

# Scale the likelihood by this factor to artificially reweight it w.r.t. the prior.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].initial.weightsMultiplier=1.0

# Minimum initial radius in pixels (used to regularize initial moments-based PSF deconvolution)
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].minInitialRadius=0.1

# Field name prefix of the Shapelet PSF approximation used to convolve the galaxy model; must contain a set of fields matching the schema defined by shapelet.MultiShapeletFunctionKey.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].psfName='modelfit_DoubleShapeletPsfApprox'

# Mask planes that indicate pixels that should be ignored in the fit.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.badMaskPlanes=['EDGE', 'SAT', 'BAD', 'NO_DATA']

# Abort if the fit region grows beyond this many pixels.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.maxArea=100000

# Maximum fraction of pixels that may be ignored due to masks; more than this and we don't even try.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.maxBadPixelFraction=0.1

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the maximum final fit region size.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMax=3.0

# Use this multiple of the initial fit ellipse then grow by the PSF width to determine the minimum final fit region size.
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nFitRadiiMin=1.0

# Use this multiple of the Kron ellipse to set the fit region (for the final fit region, subject to the nFitRadiiMin and nFitRadiiMax constraints).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nKronRadii=1.5

# Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaGrow=2.0

# If the Kron radius is less than this multiple of the PSF width, ignore it and fall back to a PSF-oriented ellipse scaled to match the area of the footprint or this radius (whichever is larger).
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].region.nPsfSigmaMin=4.0

# whether to run this plugin in single-object mode
config.processCcd.calibrate.measurement.undeblended['modelfit_CModel'].doMeasure=True

config.processCcd.calibrate.measurement.undeblended.names=[]
# Run subtask to apply aperture correction
config.processCcd.calibrate.doApCorr=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.processCcd.calibrate.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.processCcd.calibrate.applyApCorr.doFlagApCorrFailures=True

# flux measurement algorithms to be aperture-corrected by reference to another algorithm; this is a mapping alg1:alg2, where 'alg1' is the algorithm being corrected, and 'alg2' is the algorithm supplying the corrections
config.processCcd.calibrate.applyApCorr.proxies={}

# critical ratio of model to psf flux
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# correction factor for modelFlux error
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# correction factor for psfFlux error
config.processCcd.calibrate.catalogCalculation.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

config.processCcd.calibrate.catalogCalculation.plugins.names=['base_FootprintArea', 'base_ClassificationExtendedness']
# Run fake sources injection task
config.processCcd.calibrate.doInsertFakes=False

# Mask plane to set on pixels affected by fakes.  Will be added if not already present.
config.processCcd.calibrate.insertFakes.maskPlaneName='FAKE'

# Write the calexp? If fakes have been added then we do not want to write out the calexp as a normal calexp but as a fakes_calexp.
config.processCcd.calibrate.doWriteExposure=True

# name for connection icSourceSchema
config.processCcd.calibrate.connections.icSourceSchema='icSrc_schema'

# name for connection outputSchema
config.processCcd.calibrate.connections.outputSchema='src_schema'

# name for connection exposure
config.processCcd.calibrate.connections.exposure='icExp'

# name for connection background
config.processCcd.calibrate.connections.background='icExpBackground'

# name for connection icSourceCat
config.processCcd.calibrate.connections.icSourceCat='icSrc'

# name for connection astromRefCat
config.processCcd.calibrate.connections.astromRefCat='gaia_dr2_20200414'

# name for connection photoRefCat
config.processCcd.calibrate.connections.photoRefCat='ps1_pv3_3pi_20170110'

# name for connection outputExposure
config.processCcd.calibrate.connections.outputExposure='calexp'

# name for connection outputCat
config.processCcd.calibrate.connections.outputCat='src'

# name for connection outputBackground
config.processCcd.calibrate.connections.outputBackground='calexpBackground'

# name for connection matches
config.processCcd.calibrate.connections.matches='srcMatch'

# name for connection matchesDenormalized
config.processCcd.calibrate.connections.matchesDenormalized='srcMatchFull'

# Do postprocessing tasks to write parquet Source Table?
config.doMakeSourceTable=False

# ('Save the parquet version of the full src catalog?', 'Only respected if doMakeSourceTable')
config.doSaveWideSourceTable=False

# Add local photoCalib columns from the calexp.photoCalib? Should only set True if generating Source Tables from older src tables which do not already have local calib columns
config.writeSourceTable.doApplyExternalPhotoCalib=False

# Add local WCS columns from the calexp.wcs? Should only set True if generating Source Tables from older src tables which do not already have local calib columns
config.writeSourceTable.doApplyExternalSkyWcs=False

# Path to YAML file specifying functors to be computed
config.transformSourceTable.functorFile=None

# List of CCDs to ignore when processing
config.ignoreCcdList=[]

# DataId key corresponding to a single sensor
config.ccdKey='ccdnum'

