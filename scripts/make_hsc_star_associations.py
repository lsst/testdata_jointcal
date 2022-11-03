#!/usr/bin/env python
"""Script to run the isolated star association pipeline on previous hsc testdata.

This script will create a butler, run the isolated star task, and then copy the
outputs to the final location.  A temporary export file is also created which
was used as a template to add the new datasets to the existing exports.yaml.

This script only needed to be run a single time as part of DM-36726.
"""

import os
import shutil
import tempfile
import lsst.daf.butler as dafButler
import lsst.pipe.base as pipeBase
import lsst.utils
from lsst.pipe.base import Pipeline
from lsst.ctrl.mpexec import SimplePipelineExecutor
from lsst.pipe.tasks.script.registerSkymap import registerSkymap


# Set up the repo and butler.

ROOT = os.path.abspath('./')

data_dir = lsst.utils.getPackageDir("testdata_jointcal")
temp_dir = tempfile.mkdtemp(dir=ROOT, prefix="AssociationRun")
export_path = os.path.join(data_dir, "hsc/repo")
export_file = os.path.join(data_dir, "hsc", "exports.yaml")

repo = os.path.join(temp_dir, "temprepo")
_ = dafButler.Butler.makeRepo(repo)
butler = dafButler.Butler(repo, writeable=True)
instrInstance = pipeBase.Instrument.from_string("lsst.obs.subaru.HyperSuprimeCam")
instrInstance.register(butler.registry)
butler.import_(directory=export_path, filename=export_file, transfer="symlink",
               skip_dimensions={"instrument", "detector", "physical_filter"})

# Run the association pipeline.

butler = SimplePipelineExecutor.prep_butler(
    repo,
    inputs=["skymaps", "HSC/testdata"],
    output="HSC/testdata_assoc",
)

pipeline = Pipeline.fromFile(os.path.join(data_dir, "scripts", "pipelines", "isolatedStarAssociationHsc.yaml"))
executor = SimplePipelineExecutor.from_pipeline(
    pipeline,
    where="physical_filter IN ('HSC-G', 'HSC-R', 'HSC-I')",
    root=repo,
    butler=butler,
)
quanta = executor.run(register_dataset_types=True)

# Copy the output files into the right place.

starcat_file = butler.getURI('isolated_star_cat', tract=9697).ospath
starsource_file = butler.getURI('isolated_star_sources', tract=9697).ospath

out_dir = os.path.join(data_dir, "hsc", "repo", "isolated_stars")
os.makedirs(out_dir, exist_ok=True)

shutil.copy(starcat_file, out_dir)
shutil.copy(starsource_file, out_dir)

# Make a temporary export file with the dimensions

with butler.export(filename="exports_temporary_isolated.yaml") as export:
    export.saveDatasets(butler.registry.queryDatasets(collections="HSC/testdata_assoc",
                                                      datasetType="isolated_star_cat"))
    export.saveDatasets(butler.registry.queryDatasets(collections="HSC/testdata_assoc",
                                                      datasetType="isolated_star_sources"))

# By hand, update the file paths and unique ids in the export file to merge
# with the main export file. (I'm not sure how else to do this.)
