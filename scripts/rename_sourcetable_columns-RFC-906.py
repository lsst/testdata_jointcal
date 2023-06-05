#!/usr/bin/env python

"""
Script to rename source table parquet columns as part of DM-37196.
"""

import os
import glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import lsst.utils


data_dir = lsst.utils.getPackageDir("testdata_jointcal")

directories = [
    "hsc/repo/isolated_stars/",
    "hsc/repo/01291/HSC-R/output",
    "hsc/repo/01296/HSC-I/output",
    "hsc/repo/01297/HSC-I/output",
    "cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060520/r/r.MP9601/849375",
    "cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060602/r/r.MP9601/850587",
]

for directory in directories:
    parquet_files = sorted(glob.glob(os.path.join(data_dir, directory, "*.parq")))

    for parquet_file in parquet_files:
        print(f"Processing {parquet_file}")
        df = pd.read_parquet(parquet_file)

        columns = list(df.columns)
        mapper = {}
        for column in columns:
            if column == "decl":
                new_column = "dec"
                mapper[column] = new_column

        print(f"Renaming columns: {mapper}")
        df = df.rename(columns=mapper)

        os.remove(parquet_file)
        table = pa.Table.from_pandas(df)
        pq.write_table(table, parquet_file)
