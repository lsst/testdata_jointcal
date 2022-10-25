#!/usr/bin/env python

import os
import glob
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import lsst.utils


data_dir = lsst.utils.getPackageDir("testdata_jointcal")

directories = [
    'hsc/repo/01291/HSC-R/output',
    'hsc/repo/01296/HSC-I/output',
    'hsc/repo/01297/HSC-I/output',
    'cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060520/r/r.MP9601/849375',
    'cfht/repo/singleFrame/20211113T012828Z/sourceTable_visit/20060602/r/r.MP9601/850587',
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
            elif column[0].isupper():
                new_column = column[0].lower() + column[1:]
                mapper[column] = new_column

        df = df.rename(columns=mapper)

        os.remove(parquet_file)

        table = pa.Table.from_pandas(df)
        pq.write_table(table, parquet_file)
