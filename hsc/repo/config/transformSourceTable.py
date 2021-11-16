import lsst.pipe.tasks.postprocess
assert type(config)==lsst.pipe.tasks.postprocess.TransformSourceTableConfig, 'config is of type %s.%s instead of lsst.pipe.tasks.postprocess.TransformSourceTableConfig' % (type(config).__module__, type(config).__name__)
# Path to YAML file specifying functors to be computed
config.functorFile='./scripts/policy/Source_hsc.yaml'

