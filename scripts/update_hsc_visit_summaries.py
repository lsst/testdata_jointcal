#!/usr/bin/env python

import os
import glob

import lsst.afw.table


paths = [
    "hsc/repo/01297/HSC-I/output",
    "hsc/repo/01291/HSC-R/output",
    "hsc/repo/01296/HSC-I/output",
]

for path in paths:
    summfiles = sorted(glob.glob(os.path.join(path, "visitSummary*.fits")))

    for summfile in summfiles:
        summary = lsst.afw.table.ExposureCatalog.readFits(summfile)

        mapper = lsst.afw.table.SchemaMapper(summary.schema)
        mapper.addMinimalSchema(lsst.afw.table.ExposureTable.makeMinimalSchema())
        for name in summary.schema.getNames():
            if name not in mapper.getOutputSchema().getNames():
                mapper.addMapping(summary.schema[name].asKey())

        mapper.editOutputSchema().addField(
            "skyBg",
            type="F",
            doc="Average sky background (ADU)",
        )
        mapper.editOutputSchema().addField(
            "skyNoise",
            type="F",
            doc="Average sky noise (ADU)",
        )

        output_schema = mapper.getOutputSchema()

        new_summary = lsst.afw.table.ExposureCatalog(output_schema)
        new_summary.reserve(len(summary))
        new_summary.extend(summary, mapper=mapper)

        new_summary.setMetadata(summary.getMetadata())

        # We leave the skyBg as default nans because they aren't measured.
        # The important thing is to have these in the catalog.

        new_summary.writeFits(summfile)
