#!/usr/bin/env python

"""
Script to rename source table parquet columns as part of DM-31889.
This script renames gen2 ``ccd`` to modern ``detector``; renames old ``filter``
to modern ``physical_filter``; and changes columns with leading capitals (used
in gen2) to leading lower-case.
"""

import os
import glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import lsst.utils


data_dir = lsst.utils.getPackageDir("testdata_jointcal")

directories = [
    "hsc/repo/01291/HSC-R/output",
    "hsc/repo/01296/HSC-I/output",
    "hsc/repo/01297/HSC-I/output",
    "cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060520/r/r.MP9601/849375",
    "cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060602/r/r.MP9601/850587",
]

for directory in directories:
    parquet_files = sorted(glob.glob(os.path.join(data_dir, directory, "*.parq")))

    for parquet_file in parquet_files:
        df = pd.read_parquet(parquet_file)

        columns = list(df.columns)

        mapper = {}
        for column in columns:
            if column == "ccd":
                new_column = "detector"
                mapper[column] = new_column
            elif column == "filter":
                new_column = "physical_filter"
                mapper[column] = new_column
            elif column[0].isupper():
                new_column = column[0].lower() + column[1:]
                mapper[column] = new_column

        df = df.rename(columns=mapper)

        # Add band if physical_filter is in the columns.
        if "physical_filter" in columns and "band" not in columns:
            if df["physical_filter"].values[0] == "HSC-R":
                band = "r"
            elif df["physical_filter"].values[0] == "HSC-I":
                band = "i"
            elif df["physical_filter"].values[0] == "r.MP9601":
                band = "r"

            df["band"] = [band]*len(df)

        os.remove(parquet_file)

        table = pa.Table.from_pandas(df)
        pq.write_table(table, parquet_file)
