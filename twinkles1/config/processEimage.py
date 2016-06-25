import lsst.obs.lsstSim.processEimage
assert type(config)==lsst.obs.lsstSim.processEimage.ProcessEimageConfig, 'config is of type %s.%s instead of lsst.obs.lsstSim.processEimage.ProcessEimageConfig' % (type(config).__module__, type(config).__name__)
# Perform calibration?
config.doCalibrate=True

# Interpolate over defects? (ignored unless you provide a list of defects)
config.charImage.repair.doInterpolate=False

# Type of statistic to calculate edge fallbackValue for interpolation
# Allowed values:
# 	None	Field is optional
# 	MEAN	mean
# 	MEDIAN	median
# 	USER	user value set in fallbackUserValue config
# 	MEANCLIP	clipped mean
# 
config.charImage.repair.interp.fallbackValueType='MEANCLIP'

# Add a Gaussian to represent wings?
config.charImage.repair.interp.modelPsf.addWing=True

# wing amplitude, as a multiple of core amplitude (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingAmplitude=0.1

# Minimum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.minSize=5

# Maximum kernel size if using sizeFactor (pixels); ignored if size is not None
config.charImage.repair.interp.modelPsf.maxSize=None

# Kernel size as a factor of fwhm (dimensionless); size = sizeFactor * fwhm; ignored if size is not None
config.charImage.repair.interp.modelPsf.sizeFactor=3.0

# Default FWHM of Gaussian model of core of star (pixels)
config.charImage.repair.interp.modelPsf.defaultFwhm=3.0

# wing width, as a multiple of core width (dimensionless); ignored if addWing false
config.charImage.repair.interp.modelPsf.wingFwhmFactor=2.5

# Kernel size (width and height) (pixels); if None then sizeFactor is used
config.charImage.repair.interp.modelPsf.size=None

# If fallbackValueType is 'USER' then use this as the fallbackValue; ignored otherwise
config.charImage.repair.interp.fallbackUserValue=0.0

# Smoothly taper to the fallback value at the edge of the image?
config.charImage.repair.interp.useFallbackValueAtEdge=True

# Allow negative values for egde interpolation fallbackValue?  If False, set fallbackValue to max(fallbackValue, 0.0)
config.charImage.repair.interp.negativeFallbackAllowed=True

# Find and mask out cosmic rays?
config.charImage.repair.doCosmicRay=False

# Don't interpolate over CR pixels
config.charImage.repair.cosmicray.keepCRs=False

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac=2.5

# used in condition 3 for CR; see CR.cc code
config.charImage.repair.cosmicray.cond3_fac2=0.6

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.charImage.repair.cosmicray.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.charImage.repair.cosmicray.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.charImage.repair.cosmicray.background.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.charImage.repair.cosmicray.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.repair.cosmicray.background.binSize=100000

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.repair.cosmicray.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.charImage.repair.cosmicray.background.statisticsProperty='MEDIAN'

# Use Approximate (Chebyshev) to model background.
config.charImage.repair.cosmicray.background.useApprox=False

# number of times to look for contaminated pixels near known CR pixels
config.charImage.repair.cosmicray.niteration=3

# maximum number of contaminated pixels
config.charImage.repair.cosmicray.nCrPixelMax=10000

# CRs must be > this many sky-sig above sky
config.charImage.repair.cosmicray.minSigma=6.0

# CRs must have > this many DN (== electrons/gain) in initial detection
config.charImage.repair.cosmicray.min_DN=150.0

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.charImage.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.charImage.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.charImage.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.charImage.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.background.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.charImage.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.charImage.background.useApprox=False

# Write icExp and icExpBackground in addition to icSrc? Ignored if doWrite False.
config.charImage.doWriteExposure=True

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.charImage.astrometry.matcher.maxOffsetPix=300

# Type of source flux; typically one of Ap or Psf
config.charImage.astrometry.matcher.sourceFluxType='Ap'

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.charImage.astrometry.matcher.numBrightStars=50

# number of points to define a shape for matching
config.charImage.astrometry.matcher.numPointsForShape=6

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.charImage.astrometry.matcher.minSnr=40.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.charImage.astrometry.matcher.allowedNonperpDeg=3.0

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.charImage.astrometry.matcher.maxMatchDistArcSec=3.0

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.charImage.astrometry.matcher.minFracMatchedPairs=0.3

# maximum determinant of linear transformation matrix for a usable solution
config.charImage.astrometry.matcher.maxDeterminant=0.02

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.charImage.astrometry.matcher.maxRotationDeg=1.0

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.charImage.astrometry.matcher.minMatchedPairs=30

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.charImage.astrometry.minMatchDistanceArcSec=0.001

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.charImage.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.charImage.astrometry.maxIter=3

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.charImage.astrometry.matchDistanceSigma=2.0

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.charImage.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.maxScatterArcsec=10.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.numRejIter=1

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.charImage.astrometry.wcsFitter.numIter=3

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.charImage.astrometry.wcsFitter.order=4

# Number of iterations of detect sources, measure sources, estimate PSF. If useSimplePsf='all_iter' then 2 should be plenty; otherwise more may be wanted.
# 	Valid Range = [1,inf)
config.charImage.psfIterations=2

# Measure PSF? If False then keep the existing PSF model (which must exist) and use that model for all operations.
config.charImage.doMeasurePsf=True

# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	None	Field is optional
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.charImage.detectAndMeasure.measurement.doApplyApCorr='yes'

# When measuring, replace other detected footprints with noise?
config.charImage.detectAndMeasure.measurement.doReplaceWithNoise=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.charImage.detectAndMeasure.measurement.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.charImage.detectAndMeasure.measurement.applyApCorr.doFlagApCorrFailures=True

# The seed multiplier value to use for random number generation.  0 will not set seed.
config.charImage.detectAndMeasure.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.charImage.detectAndMeasure.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.charImage.detectAndMeasure.measurement.noiseReplacer.noiseSource='measure'

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_NaiveDipoleCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.charImage.detectAndMeasure.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].doMeasure=True

# Maximum flux ratio in either lobe to be considered a dipole
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].maxFluxRatio=0.65

# Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].minSn=7.0710678118654755

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.charImage.detectAndMeasure.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.charImage.detectAndMeasure.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.charImage.detectAndMeasure.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.charImage.detectAndMeasure.measurement.plugins['base_SdssCentroid'].binmax=16

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].doMeasure=True

# correction factor for psfFlux error
config.charImage.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.charImage.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.charImage.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.charImage.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.charImage.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.charImage.detectAndMeasure.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.charImage.detectAndMeasure.measurement.plugins['base_GaussianFlux'].background=0.0

# Relative weighting of pre-subtraction images (higher -> greater influence of pre-sub.
#         images on fit)
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].relWeight=0.5

# Assume dipole is not separated by more than maxSeparation * psfSigma
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxSeparation=5.0

# Set whether and how to fit for linear gradient in pre-sub. images. Possible values:
#         0: do not fit background at all
#         1 (default): pre-fit the background using linear least squares and then do not fit it as part
#             of the dipole fitting optimization
#         2: pre-fit the background using linear least squares (as in 1), and use the parameter
#             estimates from that fit as starting parameters for an integrated "re-fit" of the background
#         
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].fitBackground=1

# Maximum flux ratio in either lobe to be considered a dipole
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxFluxRatio=0.65

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].doMeasure=True

# Include parameters to fit for negative values (flux, gradient) separately from pos.
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].fitSeparateNegParams=False

# Maximum Chi2/DoF significance of fit to be considered a dipole.
#         Default value means "Choose a chi2DoF corresponding to a significance level of at most 0.05"
#         (note this is actually a significance, not a chi2 value).
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxChi2DoF=0.05

# Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].minSn=7.0710678118654755

# Fit tolerance
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].tolerance=1e-07

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_NaiveDipoleFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.charImage.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.charImage.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.charImage.detectAndMeasure.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.charImage.detectAndMeasure.measurement.plugins['base_Blendedness'].doShape=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.charImage.detectAndMeasure.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.charImage.detectAndMeasure.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.charImage.detectAndMeasure.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].doMeasure=True

# Default initial step size for flux in non-linear fitter
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].stepSizeFlux=1.0

# How many sigma the error bars of the non-linear fitter represent
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].errorDef=1.0

# Default initial step size for coordinates in non-linear fitter
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].stepSizeCoord=0.10000000149011612

# Maximum function calls for non-linear fitter; 0 = unlimited
config.charImage.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].maxFnCalls=100000

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.detectAndMeasure.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.charImage.detectAndMeasure.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.charImage.detectAndMeasure.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.charImage.detectAndMeasure.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Maximum centroid shift, limited to 2-10
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# whether to run this plugin in single-object mode
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.charImage.detectAndMeasure.measurement.plugins['base_SdssShape'].maxIter=100

config.charImage.detectAndMeasure.measurement.plugins.names=['base_CircularApertureFlux', 'base_PixelFlags', 'base_PsfFlux', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape']
# the name of the flux measurement algorithm used for calibration
config.charImage.detectAndMeasure.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.charImage.detectAndMeasure.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.charImage.detectAndMeasure.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.charImage.detectAndMeasure.measurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.charImage.detectAndMeasure.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.charImage.detectAndMeasure.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.charImage.detectAndMeasure.measurement.slots.psfFlux='base_PsfFlux'

# Deblend sources?
config.charImage.detectAndMeasure.doDeblend=False

# Estimate the background again after final source detection?
config.charImage.detectAndMeasure.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.charImage.detectAndMeasure.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.charImage.detectAndMeasure.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.charImage.detectAndMeasure.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.charImage.detectAndMeasure.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.charImage.detectAndMeasure.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.charImage.detectAndMeasure.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detectAndMeasure.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.charImage.detectAndMeasure.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.charImage.detectAndMeasure.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.charImage.detectAndMeasure.detection.includeThresholdMultiplier=10.0

# Pixels should be grown as isotropically as possible (slower)
config.charImage.detectAndMeasure.detection.isotropicGrow=False

# Fiddle factor to add to the background; debugging only
config.charImage.detectAndMeasure.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.charImage.detectAndMeasure.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.charImage.detectAndMeasure.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.charImage.detectAndMeasure.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.charImage.detectAndMeasure.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.charImage.detectAndMeasure.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.charImage.detectAndMeasure.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.charImage.detectAndMeasure.detection.background.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.charImage.detectAndMeasure.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.charImage.detectAndMeasure.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.charImage.detectAndMeasure.detection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.charImage.detectAndMeasure.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.charImage.detectAndMeasure.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.charImage.detectAndMeasure.detection.thresholdValue=5.0

# Compute aperture corrections and set ApCorrMap in exposure?
config.charImage.detectAndMeasure.doMeasureApCorr=True

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	None	Field is optional
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 
config.charImage.detectAndMeasure.deblend.edgeHandling='ramp'

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.charImage.detectAndMeasure.deblend.assignStrayFlux=True

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.charImage.detectAndMeasure.deblend.maskLimits={}

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	None	Field is optional
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 
config.charImage.detectAndMeasure.deblend.strayFluxRule='trim'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.charImage.detectAndMeasure.deblend.catchFailures=False

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	necessary	When there is not an extended object in the footprint
# 	None	Field is optional
# 
config.charImage.detectAndMeasure.deblend.strayFluxToPointSources='necessary'

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.charImage.detectAndMeasure.deblend.psfChisq2b=1.5

# Mask planes to ignore when performing statistics
config.charImage.detectAndMeasure.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.detectAndMeasure.deblend.maxFootprintArea=1000000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.detectAndMeasure.deblend.minFootprintAxisRatio=0.0

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.charImage.detectAndMeasure.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.charImage.detectAndMeasure.deblend.psfChisq1=1.5

# Find stray flux---flux not claimed by any child in the deblender.
config.charImage.detectAndMeasure.deblend.findStrayFlux=True

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.charImage.detectAndMeasure.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.charImage.detectAndMeasure.deblend.psfChisq2=1.5

# Guarantee that all peaks produce a child source.
config.charImage.detectAndMeasure.deblend.propagateAllPeaks=False

# Mask name for footprints not deblended, or None
config.charImage.detectAndMeasure.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.charImage.detectAndMeasure.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.charImage.detectAndMeasure.deblend.clipStrayFluxFraction=0.001

# Field name prefix for the flux other measurements should be aperture corrected to match
config.charImage.detectAndMeasure.measureApCorr.refFluxName='slot_CalibFlux'

# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.charImage.detectAndMeasure.measureApCorr.minDegreesOfFreedom=1

# Name of a flag field that is True for stars that should be used.
config.charImage.detectAndMeasure.measureApCorr.starSelector.field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.charImage.detectAndMeasure.measureApCorr.starSelector.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.detectAndMeasure.measureApCorr.starSelector.borderWidth=0

# size of the kernel to create
config.charImage.detectAndMeasure.measureApCorr.starSelector.kernelSize=21

# Number of standard devisations to clip at
config.charImage.detectAndMeasure.measureApCorr.numSigmaClip=3.0

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.charImage.detectAndMeasure.measureApCorr.fitConfig.triangular=True

# maximum Chebyshev function order in x
config.charImage.detectAndMeasure.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.charImage.detectAndMeasure.measureApCorr.fitConfig.orderY=2

# Number of iterations for sigma clipping
config.charImage.detectAndMeasure.measureApCorr.numIter=4

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.charImage.refObjLoader.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.charImage.refObjLoader.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.charImage.refObjLoader.filterMap={}

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.charImage.checkUnitsParseStrict='raise'

# floor for variance is lam*data
config.charImage.measurePsf.psfDeterminer['pca'].lam=0.05

# for psf candidate evaluation
config.charImage.measurePsf.psfDeterminer['pca'].reducedChi2ForPsfCandidates=3.0

# Use non-linear fitter for spatial variation of Kernel
config.charImage.measurePsf.psfDeterminer['pca'].nonLinearSpatialFit=False

# Mask blends in image?
config.charImage.measurePsf.psfDeterminer['pca'].doMaskBlends=True

# size of cell used to determine PSF (pixels, column direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellX=256

# size of cell used to determine PSF (pixels, row direction)
config.charImage.measurePsf.psfDeterminer['pca'].sizeCellY=256

# Rejection threshold (stdev) for candidates based on spatial fit
config.charImage.measurePsf.psfDeterminer['pca'].spatialReject=2.0

# radius of the kernel to create, relative to the square root of the stellar quadrupole moments
config.charImage.measurePsf.psfDeterminer['pca'].kernelSize=10.0

# Reject candidates that are blended?
config.charImage.measurePsf.psfDeterminer['pca'].doRejectBlends=False

# number of iterations of PSF candidate star list
config.charImage.measurePsf.psfDeterminer['pca'].nIterForPsf=0

# Should each PSF candidate be given the same weight, independent of magnitude?
config.charImage.measurePsf.psfDeterminer['pca'].constantWeight=True

# Maximum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMax=45

# number of stars per psf Cell for spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCellSpatialFit=5

# Number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.psfDeterminer['pca'].borderWidth=0

# number of eigen components for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nEigenComponents=4

# Threshold (stdev) for rejecting extraneous pixels around candidate; applied if positive
config.charImage.measurePsf.psfDeterminer['pca'].pixelThreshold=0.0

# specify spatial order for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].spatialOrder=2

# tolerance of spatial fitting
config.charImage.measurePsf.psfDeterminer['pca'].tolerance=0.01

# Minimum radius of the kernel
config.charImage.measurePsf.psfDeterminer['pca'].kernelSizeMin=25

# number of stars per psf cell for PSF kernel creation
config.charImage.measurePsf.psfDeterminer['pca'].nStarPerCell=3

config.charImage.measurePsf.psfDeterminer.name='pca'
# Standard deviation of width allowed to be interpreted as good stars
config.charImage.measurePsf.starSelector.widthStdAllowed=1.0

# maximum width to include in histogram
config.charImage.measurePsf.starSelector.widthMax=10.0

# minimum width to include in histogram
config.charImage.measurePsf.starSelector.widthMin=0.0

# size of the kernel to create
config.charImage.measurePsf.starSelector.kernelSize=21

# specify the minimum psfFlux for good Psf Candidates
config.charImage.measurePsf.starSelector.fluxMin=12500.0

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.charImage.measurePsf.starSelector.borderWidth=0

# List of flags which cause a source to be rejected as bad
config.charImage.measurePsf.starSelector.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# specify the maximum psfFlux for good Psf Candidates (ignored if == 0)
config.charImage.measurePsf.starSelector.fluxMax=0.0

# Name of field in Source to use for flux measurement
config.charImage.measurePsf.starSelector.sourceFluxField='base_GaussianFlux_flux'

# Keep objects within this many sigma of cluster 0's median
config.charImage.measurePsf.starSelector.nSigmaClip=2.0

# This number will be multplied by the exposure ID to set the random seed for reserving candidates
config.charImage.measurePsf.reserveSeed=1

# Fraction of PSF candidates to reserve from fitting; none if <= 0
config.charImage.measurePsf.reserveFraction=-1.0

# Width and height of PSF model, in pixels. Must be odd.
# 	Valid Range = [1,inf)
config.charImage.installSimplePsf.width=11

# Estimated FWHM of simple Gaussian PSF model, in pixels. Ignored if input exposure has a PSF model.
config.charImage.installSimplePsf.fwhm=3.5322300675464238

# Persist results?
config.charImage.doWrite=True

# Replace the existing PSF model with a simplified version that has the same sigma at the start of each PSF determination iteration? Doing so makes PSF determination converge more robustly and quickly.
config.charImage.useSimplePsf=True

# Seed for random number generator
config.rngSeed=1234567890

# Perform astrometric calibration?
config.calibrate.doAstrometry=True

# Include HeavyFootprint data in source table? If false then heavy footprints are saved as normal footprints, which saves some space
config.calibrate.doWriteHeavyFootprintsInSources=True

# Write reference matches (ignored if doWrite false)?
config.calibrate.doWriteMatches=True

# Apply aperture corrections? Silently ignored if endOrder <= lsst.meas.base.APCORR_ORDER when calling run
# Allowed values:
# 	yes	apply aperture corrections; fail if data not available
# 	noButWarn	do not apply aperture corrections, but warn if data available (since aperture corrections could have been applied)
# 	None	Field is optional
# 	yesOrWarn	apply aperture corrections if data available, else warn
# 	no	do not apply aperture corrections
# 
config.calibrate.detectAndMeasure.measurement.doApplyApCorr='yes'

# When measuring, replace other detected footprints with noise?
config.calibrate.detectAndMeasure.measurement.doReplaceWithNoise=True

# flux measurement algorithms in getApCorrNameSet() to ignore; if a name is listed that does not appear in getApCorrNameSet() then a warning is logged
config.calibrate.detectAndMeasure.measurement.applyApCorr.ignoreList=[]

# set the general failure flag for a flux when it cannot be aperture-corrected?
config.calibrate.detectAndMeasure.measurement.applyApCorr.doFlagApCorrFailures=True

# The seed multiplier value to use for random number generation.  0 will not set seed.
config.calibrate.detectAndMeasure.measurement.noiseReplacer.noiseSeedMultiplier=1

# Add ann offset to the generated noise.
config.calibrate.detectAndMeasure.measurement.noiseReplacer.noiseOffset=0.0

# How to choose mean and variance of the Gaussian noise we generate?
# Allowed values:
# 	variance	Mean = 0, variance = the image's variance
# 	meta	Mean = 0, variance = the "BGMEAN" metadata entry
# 	measure	Measure clipped mean and variance from the whole image
# 
config.calibrate.detectAndMeasure.measurement.noiseReplacer.noiseSource='measure'

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_NaiveDipoleCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_FPPosition'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_NaiveCentroid'].doMeasure=True

# Value to subtract from the image pixel values
config.calibrate.detectAndMeasure.measurement.plugins['base_NaiveCentroid'].background=0.0

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_InputCount'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].doMeasure=True

# Maximum flux ratio in either lobe to be considered a dipole
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].maxFluxRatio=0.65

# Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_ClassificationDipole'].minSn=7.0710678118654755

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_PeakLikelihoodFlux'].doMeasure=True

# Name of warping kernel (e.g. "lanczos4") used to compute the peak
config.calibrate.detectAndMeasure.measurement.plugins['base_PeakLikelihoodFlux'].warpingKernelName='lanczos4'

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssCentroid'].doMeasure=True

# if the peak's less than this insist on binning at least once
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssCentroid'].peakMin=-1.0

# fiddle factor for adjusting the binning
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssCentroid'].wfac=1.5

# maximum allowed binning
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssCentroid'].binmax=16

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].doMeasure=True

# correction factor for psfFlux error
config.calibrate.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].psfErrFactor=0.0

# correction factor for modelFlux error
config.calibrate.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].modelErrFactor=0.0

# critical ratio of model to psf flux
config.calibrate.detectAndMeasure.measurement.plugins['base_ClassificationExtendedness'].fluxRatio=0.925

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].shiftKernel='lanczos5'

# Radius (in pixels) of apertures.
config.calibrate.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0]

# Maximum radius (in pixels) for which the sinc algorithm should be used instead of the faster naive algorithm.  For elliptical apertures, this is the minor axis radius.
config.calibrate.detectAndMeasure.measurement.plugins['base_CircularApertureFlux'].maxSincRadius=10.0

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_Jacobian'].doMeasure=True

# Nominal pixel size (arcsec)
config.calibrate.detectAndMeasure.measurement.plugins['base_Jacobian'].pixelScale=0.5

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_GaussianCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_GaussianFlux'].doMeasure=True

# FIXME! NEVER DOCUMENTED!
config.calibrate.detectAndMeasure.measurement.plugins['base_GaussianFlux'].background=0.0

# Relative weighting of pre-subtraction images (higher -> greater influence of pre-sub.
#         images on fit)
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].relWeight=0.5

# Assume dipole is not separated by more than maxSeparation * psfSigma
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxSeparation=5.0

# Set whether and how to fit for linear gradient in pre-sub. images. Possible values:
#         0: do not fit background at all
#         1 (default): pre-fit the background using linear least squares and then do not fit it as part
#             of the dipole fitting optimization
#         2: pre-fit the background using linear least squares (as in 1), and use the parameter
#             estimates from that fit as starting parameters for an integrated "re-fit" of the background
#         
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].fitBackground=1

# Maximum flux ratio in either lobe to be considered a dipole
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxFluxRatio=0.65

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].doMeasure=True

# Include parameters to fit for negative values (flux, gradient) separately from pos.
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].fitSeparateNegParams=False

# Maximum Chi2/DoF significance of fit to be considered a dipole.
#         Default value means "Choose a chi2DoF corresponding to a significance level of at most 0.05"
#         (note this is actually a significance, not a chi2 value).
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].maxChi2DoF=0.05

# Minimum quadrature sum of positive+negative lobe S/N to be considered a dipole
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].minSn=7.0710678118654755

# Fit tolerance
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_DipoleFit'].tolerance=1e-07

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_NaiveDipoleFlux'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].doMeasure=True

# Warping kernel used to shift Sinc photometry coefficients to different center positions
config.calibrate.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].shiftKernel='lanczos5'

# Scaling factor of PSF FWHM for aperture radius.
config.calibrate.detectAndMeasure.measurement.plugins['base_ScaledApertureFlux'].scale=3.14

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_PsfFlux'].doMeasure=True

# Mask planes that indicate pixels that should be excluded from the fit
config.calibrate.detectAndMeasure.measurement.plugins['base_PsfFlux'].badMaskPlanes=[]

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_Blendedness'].doMeasure=True

# Whether to compute quantities related to the Gaussian-weighted shape
config.calibrate.detectAndMeasure.measurement.plugins['base_Blendedness'].doShape=True

# Whether to compute HeavyFootprint dot products (the old deblend.blendedness parameter)
config.calibrate.detectAndMeasure.measurement.plugins['base_Blendedness'].doOld=True

# Whether to compute quantities related to the Gaussian-weighted flux
config.calibrate.detectAndMeasure.measurement.plugins['base_Blendedness'].doFlux=True

# Radius factor that sets the maximum extent of the weight function (and hence the flux measurements)
config.calibrate.detectAndMeasure.measurement.plugins['base_Blendedness'].nSigmaWeightMax=3.0

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_PeakCentroid'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].doMeasure=True

# Default initial step size for flux in non-linear fitter
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].stepSizeFlux=1.0

# How many sigma the error bars of the non-linear fitter represent
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].errorDef=1.0

# Default initial step size for coordinates in non-linear fitter
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].stepSizeCoord=0.10000000149011612

# Maximum function calls for non-linear fitter; 0 = unlimited
config.calibrate.detectAndMeasure.measurement.plugins['ip_diffim_PsfDipoleFlux'].maxFnCalls=100000

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_PixelFlags'].doMeasure=True

# List of mask planes to be searched for which occur in the center of a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.detectAndMeasure.measurement.plugins['base_PixelFlags'].masksFpCenter=[]

# List of mask planes to be searched for which occur anywhere within a footprint. If any of the planes are found they will have a corresponding pixel flag set.
config.calibrate.detectAndMeasure.measurement.plugins['base_PixelFlags'].masksFpAnywhere=[]

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_SkyCoord'].doMeasure=True

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_Variance'].doMeasure=True

# Scale factor to apply to shape for aperture
config.calibrate.detectAndMeasure.measurement.plugins['base_Variance'].scale=5.0

# Mask planes to ignore
config.calibrate.detectAndMeasure.measurement.plugins['base_Variance'].mask=['DETECTED', 'DETECTED_NEGATIVE', 'BAD', 'SAT']

# Maximum centroid shift, limited to 2-10
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].maxShift=0.0

# Convergence tolerance for FWHM
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].tol2=9.999999747378752e-05

# Convergence tolerance for e1,e2
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].tol1=9.999999747378752e-06

# whether to run this plugin in single-object mode
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].doMeasure=True

# Additional value to add to background
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].background=0.0

# Maximum number of iterations
config.calibrate.detectAndMeasure.measurement.plugins['base_SdssShape'].maxIter=100

config.calibrate.detectAndMeasure.measurement.plugins.names=['base_CircularApertureFlux', 'base_NaiveCentroid', 'base_PixelFlags', 'base_SkyCoord', 'base_PsfFlux', 'base_Variance', 'base_GaussianCentroid', 'base_SdssCentroid', 'base_GaussianFlux', 'base_SdssShape', 'base_ClassificationExtendedness']
# the name of the flux measurement algorithm used for calibration
config.calibrate.detectAndMeasure.measurement.slots.calibFlux='base_CircularApertureFlux_12_0'

# the name of the algorithm used to set the source aperture flux slot
config.calibrate.detectAndMeasure.measurement.slots.apFlux='base_CircularApertureFlux_3_0'

# the name of the algorithm used to set the source inst flux slot
config.calibrate.detectAndMeasure.measurement.slots.instFlux='base_GaussianFlux'

# the name of the algorithm used to set source moments parameters
config.calibrate.detectAndMeasure.measurement.slots.shape='base_SdssShape'

# the name of the centroiding algorithm used to set source x,y
config.calibrate.detectAndMeasure.measurement.slots.centroid='base_SdssCentroid'

# the name of the algorithm used to set the source model flux slot
config.calibrate.detectAndMeasure.measurement.slots.modelFlux='base_GaussianFlux'

# the name of the algorithm used to set the source psf flux slot
config.calibrate.detectAndMeasure.measurement.slots.psfFlux='base_PsfFlux'

# Deblend sources?
config.calibrate.detectAndMeasure.doDeblend=True

# Estimate the background again after final source detection?
config.calibrate.detectAndMeasure.detection.reEstimateBackground=True

# Grow detections by nSigmaToGrow * sigma; if 0 then do not grow
config.calibrate.detectAndMeasure.detection.nSigmaToGrow=2.4

# detected sources with fewer than the specified number of pixels will be ignored
# 	Valid Range = [0,inf)
config.calibrate.detectAndMeasure.detection.minPixels=1

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.tempLocalBackground.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.tempLocalBackground.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.detectAndMeasure.detection.tempLocalBackground.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.detectAndMeasure.detection.tempLocalBackground.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.detectAndMeasure.detection.tempLocalBackground.algorithm='AKIMA_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.detectAndMeasure.detection.tempLocalBackground.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detectAndMeasure.detection.tempLocalBackground.binSize=64

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.tempLocalBackground.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.detectAndMeasure.detection.tempLocalBackground.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.calibrate.detectAndMeasure.detection.tempLocalBackground.useApprox=False

# Include threshold relative to thresholdValue
# 	Valid Range = [0.0,inf)
config.calibrate.detectAndMeasure.detection.includeThresholdMultiplier=1.0

# Pixels should be grown as isotropically as possible (slower)
config.calibrate.detectAndMeasure.detection.isotropicGrow=False

# Fiddle factor to add to the background; debugging only
config.calibrate.detectAndMeasure.detection.adjustBackground=0.0

# Do temporary interpolated background subtraction before footprint detection?
config.calibrate.detectAndMeasure.detection.doTempLocalBackground=False

# specifies whether to detect positive, or negative sources, or both
# Allowed values:
# 	positive	detect only positive sources
# 	negative	detect only negative sources
# 	both	detect both positive and negative sources
# 
config.calibrate.detectAndMeasure.detection.thresholdPolarity='positive'

# Approximation order in X for background Chebyshev (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.background.approxOrderX=6

# Approximation order in Y for background Chebyshev (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.background.approxOrderY=-1

# Names of mask planes to ignore while estimating the background
config.calibrate.detectAndMeasure.detection.background.ignoredPixelMask=['BAD', 'EDGE', 'DETECTED', 'DETECTED_NEGATIVE', 'NO_DATA']

# behaviour if there are too few points in grid for requested interpolation style
# Allowed values:
# 	INCREASE_NXNYSAMPLE	Increase the number of samples used to make the interpolation grid.
# 	None	Field is optional
# 	THROW_EXCEPTION	throw an exception if there are too few points
# 	REDUCE_INTERP_ORDER	use an interpolation style with a lower order.
# 
config.calibrate.detectAndMeasure.detection.background.undersampleStyle='REDUCE_INTERP_ORDER'

# how to interpolate the background values. This maps to an enum; see afw::math::Background
# Allowed values:
# 	NONE	No background estimation is to be attempted
# 	None	Field is optional
# 	LINEAR	Use linear interpolation
# 	AKIMA_SPLINE	higher-level nonlinear spline that is more robust to outliers
# 	CONSTANT	Use a single constant value
# 	NATURAL_SPLINE	cubic spline with zero second derivative at endpoints
# 
config.calibrate.detectAndMeasure.detection.background.algorithm='NATURAL_SPLINE'

# Ignore NaNs when estimating the background
config.calibrate.detectAndMeasure.detection.background.isNanSafe=False

# how large a region of the sky should be used for each background point
# 	Valid Range = [1,inf)
config.calibrate.detectAndMeasure.detection.background.binSize=256

# Use inverse variance weighting in calculation (valid only with useApprox=True)
config.calibrate.detectAndMeasure.detection.background.weighting=True

# type of statistic to use for grid points
# Allowed values:
# 	None	Field is optional
# 	MEDIAN	median
# 	MEANCLIP	clipped mean
# 	MEAN	unclipped mean
# 
config.calibrate.detectAndMeasure.detection.background.statisticsProperty='MEANCLIP'

# Use Approximate (Chebyshev) to model background.
config.calibrate.detectAndMeasure.detection.background.useApprox=False

# Grow detections to set the image mask bits, but return the original (not-grown) footprints
config.calibrate.detectAndMeasure.detection.returnOriginalFootprints=False

# specifies the desired flavor of Threshold
# Allowed values:
# 	pixel_stdev	threshold applied to per-pixel std deviation
# 	variance	threshold applied to image variance
# 	value	threshold applied to image value
# 	stdev	threshold applied to image std deviation
# 
config.calibrate.detectAndMeasure.detection.thresholdType='stdev'

# Threshold for footprints
# 	Valid Range = [0.0,inf)
config.calibrate.detectAndMeasure.detection.thresholdValue=5.0

# Compute aperture corrections and set ApCorrMap in exposure?
config.calibrate.detectAndMeasure.doMeasureApCorr=False

# What to do when a peak to be deblended is close to the edge of the image
# Allowed values:
# 	None	Field is optional
# 	ramp	Ramp down flux at the image edge by the PSF
# 	noclip	Ignore the edge when building the symmetric template.
# 	clip	Clip the template at the edge AND the mirror of the edge.
# 
config.calibrate.detectAndMeasure.deblend.edgeHandling='ramp'

# Assign stray flux to deblend children.  Implies findStrayFlux.
config.calibrate.detectAndMeasure.deblend.assignStrayFlux=True

# Mask planes with the corresponding limit on the fraction of masked pixels. Sources violating this limit will not be deblended.
config.calibrate.detectAndMeasure.deblend.maskLimits={}

# How to split flux among peaks
# Allowed values:
# 	trim	Shrink the parent footprint to pixels that are not assigned to children
# 	r-to-peak	~ 1/(1+R^2) to the peak
# 	None	Field is optional
# 	r-to-footprint	~ 1/(1+R^2) to the closest pixel in the footprint.  CAUTION: this can be computationally expensive on large footprints!
# 	nearest-footprint	Assign 100% to the nearest footprint (using L-1 norm aka Manhattan distance)
# 
config.calibrate.detectAndMeasure.deblend.strayFluxRule='trim'

# If True, catch exceptions thrown by the deblender, log them, and set a flag on the parent, instead of letting them propagate up
config.calibrate.detectAndMeasure.deblend.catchFailures=False

# When the deblender should attribute stray flux to point sources
# Allowed values:
# 	always	Always
# 	never	Never; stray flux will not be attributed to any deblended child if the deblender thinks all peaks look like point sources
# 	necessary	When there is not an extended object in the footprint
# 	None	Field is optional
# 
config.calibrate.detectAndMeasure.deblend.strayFluxToPointSources='necessary'

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (shifted PSF model #2)
config.calibrate.detectAndMeasure.deblend.psfChisq2b=1.5

# Mask planes to ignore when performing statistics
config.calibrate.detectAndMeasure.deblend.maskPlanes=['SAT', 'INTRP', 'NO_DATA']

# Maximum area for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.detectAndMeasure.deblend.maxFootprintArea=1000000

# Minimum axis ratio for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.detectAndMeasure.deblend.minFootprintAxisRatio=0.0

# Maximum linear dimension for footprints before they are ignored as large; non-positive means no threshold applied
config.calibrate.detectAndMeasure.deblend.maxFootprintSize=0

# Chi-squared per DOF cut for deciding a source is a PSF during deblending (un-shifted PSF model)
config.calibrate.detectAndMeasure.deblend.psfChisq1=1.5

# Find stray flux---flux not claimed by any child in the deblender.
config.calibrate.detectAndMeasure.deblend.findStrayFlux=True

# Footprints smaller in width or height than this value will be ignored; minimum of 2 due to PSF gradient calculation.
# 	Valid Range = [2,inf)
config.calibrate.detectAndMeasure.deblend.tinyFootprintSize=2

# Chi-squared per DOF cut for deciding a source is PSF during deblending (shifted PSF model)
config.calibrate.detectAndMeasure.deblend.psfChisq2=1.5

# Guarantee that all peaks produce a child source.
config.calibrate.detectAndMeasure.deblend.propagateAllPeaks=False

# Mask name for footprints not deblended, or None
config.calibrate.detectAndMeasure.deblend.notDeblendedMask='NOT_DEBLENDED'

# Only deblend the brightest maxNumberOfPeaks peaks in the parent (<= 0: unlimited)
config.calibrate.detectAndMeasure.deblend.maxNumberOfPeaks=0

# When splitting stray flux, clip fractions below this value to zero.
config.calibrate.detectAndMeasure.deblend.clipStrayFluxFraction=0.001

# Field name prefix for the flux other measurements should be aperture corrected to match
config.calibrate.detectAndMeasure.measureApCorr.refFluxName='slot_CalibFlux'

# Minimum number of degrees of freedom (# of valid data points - # of parameters); if this is exceeded, the order of the fit is decreased (in both dimensions), and if we can't decrease it enough, we'll raise ValueError.
# 	Valid Range = [1,inf)
config.calibrate.detectAndMeasure.measureApCorr.minDegreesOfFreedom=1

# Name of a flag field that is True for stars that should be used.
config.calibrate.detectAndMeasure.measureApCorr.starSelector.field='calib_psfUsed'

# List of flags which cause a source to be rejected as bad
config.calibrate.detectAndMeasure.measureApCorr.starSelector.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolatedCenter', 'base_PixelFlags_flag_saturatedCenter', 'base_PixelFlags_flag_crCenter', 'base_PixelFlags_flag_bad', 'base_PixelFlags_flag_interpolated']

# number of pixels to ignore around the edge of PSF candidate postage stamps
config.calibrate.detectAndMeasure.measureApCorr.starSelector.borderWidth=0

# size of the kernel to create
config.calibrate.detectAndMeasure.measureApCorr.starSelector.kernelSize=21

# Number of standard devisations to clip at
config.calibrate.detectAndMeasure.measureApCorr.numSigmaClip=3.0

# if true, only include terms where the sum of the x and y order is less than or equal to max(orderX, orderY)
config.calibrate.detectAndMeasure.measureApCorr.fitConfig.triangular=True

# maximum Chebyshev function order in x
config.calibrate.detectAndMeasure.measureApCorr.fitConfig.orderX=2

# maximum Chebyshev function order in y
config.calibrate.detectAndMeasure.measureApCorr.fitConfig.orderY=2

# Number of iterations for sigma clipping
config.calibrate.detectAndMeasure.measureApCorr.numIter=4

# Maximum allowed shift of WCS, due to matching (pixel)
# 	Valid Range = [-inf,4000)
config.calibrate.astrometry.matcher.maxOffsetPix=300

# Type of source flux; typically one of Ap or Psf
config.calibrate.astrometry.matcher.sourceFluxType='Ap'

# Number of bright stars to use
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.numBrightStars=50

# number of points to define a shape for matching
config.calibrate.astrometry.matcher.numPointsForShape=6

# Minimum allowed signal-to-noise ratio for sources used for matching (in the flux specified by sourceFluxType); <=0 for no limit
config.calibrate.astrometry.matcher.minSnr=40.0

# Allowed non-perpendicularity of x and y (degree)
# 	Valid Range = [-inf,45.0)
config.calibrate.astrometry.matcher.allowedNonperpDeg=3.0

# Maximum separation between reference objects and sources beyond which they will not be considered a match (arcsec)
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matcher.maxMatchDistArcSec=3.0

# Minimum number of matched pairs as a fraction of the smaller of the number of reference stars or the number of good sources; the actual minimum is the smaller of this value or minMatchedPairs
# 	Valid Range = [0,1)
config.calibrate.astrometry.matcher.minFracMatchedPairs=0.3

# maximum determinant of linear transformation matrix for a usable solution
config.calibrate.astrometry.matcher.maxDeterminant=0.02

# Rotation angle allowed between sources and position reference objects (degrees)
# 	Valid Range = [-inf,6.0)
config.calibrate.astrometry.matcher.maxRotationDeg=1.0

# Minimum number of matched pairs; see also minFracMatchedPairs
# 	Valid Range = [2,inf)
config.calibrate.astrometry.matcher.minMatchedPairs=30

# the match distance below which further iteration is pointless (arcsec); ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.minMatchDistanceArcSec=0.001

# If True then load reference objects and match sources but do not fit a WCS;  this simply controls whether 'run' calls 'solve' or 'loadAndMatch'
config.calibrate.astrometry.forceKnownWcs=False

# maximum number of iterations of match sources and fit WCSignored if not fitting a WCS
# 	Valid Range = [1,inf)
config.calibrate.astrometry.maxIter=3

# the maximum match distance is set to  mean_match_distance + matchDistanceSigma*std_dev_match_distance; ignored if not fitting a WCS
# 	Valid Range = [0,inf)
config.calibrate.astrometry.matchDistanceSigma=2.0

# Number of standard deviations for clipping level
# 	Valid Range = [0.0,inf)
config.calibrate.astrometry.wcsFitter.rejSigma=3.0

# maximum median scatter of a WCS fit beyond which the fit fails (arcsec); be generous, as this is only intended to catch catastrophic failures
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.maxScatterArcsec=10.0

# number of rejection iterations
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.numRejIter=1

# number of iterations of fitter (which fits X and Y separately, and so benefits from a few iterations
# 	Valid Range = [1,inf)
config.calibrate.astrometry.wcsFitter.numIter=3

# order of SIP polynomial
# 	Valid Range = [0,inf)
config.calibrate.astrometry.wcsFitter.order=4

# Fields to copy from the icSource catalog to the output catalog for matching sources Any missing fields will trigger a RuntimeError exception. If detectAndMeasure.doMeasureApCorr is True and detectAndMeasure cannot determine its own suitable candidates, then this list must include config.detectAndMeasure.measureApCorr.inputFilterFlag. Ignored if icSourceCat is not provided.
config.calibrate.icSourceFieldsToCopy=['calib_psfCandidate', 'calib_psfUsed', 'calib_psfReserved']

# Perform phometric calibration?
config.calibrate.doPhotoCal=True

# Default reference catalog filter to use if filter not specified in exposure; if blank then filter must be specified in exposure
config.calibrate.refObjLoader.defaultFilter=''

# Padding to add to 4 all edges of the bounding box (pixels)
# 	Valid Range = [0,inf)
config.calibrate.refObjLoader.pixelMargin=50

# Mapping of camera filter name: reference catalog filter name; each reference filter must exist
config.calibrate.refObjLoader.filterMap={}

# Strictness of Astropy unit compatibility check, can be 'raise', 'warn' or 'silent'
config.calibrate.checkUnitsParseStrict='raise'

# Raise an exception if photoCal fails? Ignored if doPhotoCal false.
config.calibrate.requirePhotoCal=True

# Save calibration results?
config.calibrate.doWrite=True

# Use the extendedness parameter to select objects to use in photometric calibration?
# This applies only to the sources detected on the exposure, not the reference catalog
config.calibrate.photoCal.doSelectUnresolved=False

# number of iterations
config.calibrate.photoCal.nIter=20

# Name of the source flux field to use.  The associated flag field
# ('<name>_flags') will be implicitly included in badFlags.
config.calibrate.photoCal.fluxField='slot_CalibFlux_flux'

config.calibrate.photoCal.colorterms.data={}
# Name of photometric reference catalog; used to select a color term dict in colorterms. see also applyColorTerms
config.calibrate.photoCal.photoCatName=None

# Don't use objects fainter than this magnitude
config.calibrate.photoCal.magLimit=22.0

# maximum sigma to use when clipping
config.calibrate.photoCal.sigmaMax=0.25

# List of source flag fields that must be set for a source to be used.
config.calibrate.photoCal.goodFlags=[]

# Additional magnitude uncertainty to be added in quadrature with measurement errors.
# 	Valid Range = [0.0,inf)
config.calibrate.photoCal.magErrFloor=0.0

# clip at nSigma
config.calibrate.photoCal.nSigma=3.0

# List of source flag fields that will cause a source to be rejected when they are set.
config.calibrate.photoCal.badFlags=['base_PixelFlags_flag_edge', 'base_PixelFlags_flag_interpolated', 'base_PixelFlags_flag_saturated']

# use median instead of mean to compute zeropoint
config.calibrate.photoCal.useMedian=True

# Write a field name astrom_usedByPhotoCal to the schema
config.calibrate.photoCal.doWriteOutput=True

# Apply photometric color terms to reference stars? One of:
# None: apply if colorterms and photoCatName are not None;
#       fail if color term data is not available for the specified ref catalog and filter.
# True: always apply colorterms; fail if color term data is not available for the
#       specified reference catalog and filter.
# False: do not apply.
config.calibrate.photoCal.applyColorTerms=None

# Match radius for matching icSourceCat objects to sourceCat objects (pixels)
config.calibrate.matchRadiusPix=3.0

# Raise an exception if astrometry fails? Ignored if doAstrometry false.
config.calibrate.requireAstrometry=True

import lsst.obs.lsstSim.eimageIsr
config.isr.retarget(target=lsst.obs.lsstSim.eimageIsr.EimageIsrTask, ConfigClass=lsst.obs.lsstSim.eimageIsr.EimageIsrConfig)
# Random number seed used when adding noise (passed directly to numpy at task initialization)
config.isr.rngSeed=None

# Value at which to detect saturation
config.isr.sat_val=100000

# Set mask to EDGE for a border of x pixels
config.isr.maskEdgeBorder=0

# Size of interpolation kernel in arcsec
config.isr.interp_size=0.5

# Mean of the Poisson distribution in counts
config.isr.noiseValue=1000

# Set the variance plane in the eimage?
config.isr.doSetVariance=True

# Choose method for setting the variance
# Allowed values:
# 	image	set variance from image plane
# 	None	Field is optional
# 	value	set variance to a value
# 
config.isr.varianceType='image'

# Add a flat Poisson noise background to the eimage?
config.isr.doAddNoise=False

# Value to use in the variance plane.
config.isr.varianceValue=0.01

