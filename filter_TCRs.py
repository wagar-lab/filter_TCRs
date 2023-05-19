# filter_TCRs.py
# This is the program our lab uses to QC-filter the TCRs from our input. 
# We assume that the input is a 10x file of a single donor from both compartments that has undergone single-cell sequencing filtering and has had the GEX annotations integrated.
# Requires NumPy and Pandas (use the latest version)

import sys
import time
import argparse
import pandas as pd
from numpy import savetxt

def log(message, stream=sys.stderr):
    print(f"[{time.strftime('%H:%M:%S', time.localtime())}] {message}", file=stream)


def filter_TCRs(df):
    '''
    Removes cells from the df that do not have an alpha or beta chain, or have more than one beta chain.
    '''
    # Remove rows that do not have enough data
    df = df.dropna()

    # Remove cells that do not have an alpha or beta chain
    TRA_df = df.loc[df["chain"] == "TRA"]
    TRB_df = df.loc[df["chain"] == "TRB"]
    paired = TRB_df.loc[TRB_df.barcode.isin(TRA_df.barcode)]
    df = df[df.barcode.isin(paired.barcode)]

    # Remove cells with more than one beta chain
    multiple_TRBs_df = TRB_df[TRB_df.duplicated(subset=["barcode", "chain"], keep=False)]
    df = df[~df.barcode.isin(multiple_TRBs_df)]
    return df

if __name__ == "__main__":
    # 0) Argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help=".csv file to read from", required=True)
    parser.add_argument("-o", "--outfile", type=str, help="file to save to", required=True)
    parser.add_argument("-n", "--sample", type=int, help="sample N cells from the input file")    
    args = parser.parse_args()

    # 1) Read the file.
    log(f"Reading {args.file}")
    df = pd.read_csv(args.file, engine="pyarrow")
    cells = pd.Series(df.barcode.unique())
    log(f"There are {len(df)} rows and {len(cells)} cells")
    if args.sample:
        if len(cells) < args.sample:
            log(f"Sample size bigger than cell count, sampling n={len(cells)} instead")
            sample = cells.sample(len(cells))
        else:
            log(f"Sampling {sample_size} cells from the data")
            sample = cells.sample(args.sample)
        df = df[df.barcode.isin(sample)]
    
    # 2) Split into PBMC and Tonsil (avoids problems where a cell may be high-quality in one compartment, but not the other).
    pbmc_df = df.loc[df.compartment=="pbmc"]
    tonsil_df = df.loc[df.compartment=="tonsil"]

    # 3) Filter the data and merge them together.
    pbmc_df = filter_TCRs(pbmc_df)
    tonsil_df = filter_TCRs(tonsil_df)
    df = pd.concat([pbmc_df, tonsil_df], ignore_index=True)
    filtered_cells = pd.Series(df.barcode.unique())
    log(f"There are {len(df)} rows and {len(filtered_cells)} cells")

    # 4) Write it back out
    log(f"Writing to {args.outfile}")
    header_fmt = "%s,%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%s,%s,%d,%d,%s,%s"
    savetxt(args.outfile, df.values, fmt=header_fmt, header=','.join(df.columns), comments='')
    log(f"Done writing to {args.outfile}")

