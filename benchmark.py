# benchmark.py
# Code used to produce the asymptotic performance graph from my poster.
# The input is a 10x file that has gone through single-cell sequencing filtering and has had the GEX annotations integrated.
# Unfortunately, I am not allowed to share the input data.
#
# Please read filter_TCRs.py for more details about the program.

from filter_TCRs import log, filter_TCRs
import argparse
import time
import pandas as pd

 
if __name__ == "__main__":
    # Argument parser init
    parser = argparse.ArgumentParser(description="Asymptotic performance benchmark. Measures the elapsed time of the filtering step at varying input sizes (powers of two).")
    parser.add_argument("-f", "--file", type=str, help=".csv file to read from", required=True)
    parser.add_argument("-o", "--outfile", type=str, help="file to save to")
    args = parser.parse_args()
    
    # 1) Read the file
    log(f"Reading {args.file}")
    df = pd.read_csv(args.file, engine='pyarrow')
    cells = pd.Series(df.barcode.unique())
    log(f"There are {len(df)} rows and {len(cells)} cells")
    if args.outfile:
        f = open(args.outfile, "w")
    # 2) Benchmark
    times = []
    sample_size = 2**6
    while sample_size < len(cells):
        log(f"Sampling {sample_size} cells from the data")
        sample = cells.sample(sample_size)
        df_subset = df[df.barcode.isin(sample)] 
        
        start = time.process_time_ns()
        
        # Split by donor, and then filter blood and tonsil separately
        updated_dfs = []
        for donor_name, donor_df in df_subset.groupby("donor"): 
            pbmc_df = donor_df.loc[donor_df["compartment"]=="pbmc"]
            tonsil_df = donor_df.loc[donor_df["compartment"]=="tonsil"]

            pbmc_df = filter_TCRs(pbmc_df)
            tonsil_df = filter_TCRs(tonsil_df)
            updated_df = pd.concat([pbmc_df, tonsil_df], ignore_index=True)
            updated_dfs.append(updated_df)
        final_df = pd.concat(updated_dfs).reset_index(drop=True) # Final, filtered df
        end = time.process_time_ns()
        elapsed = end - start
        times.append(elapsed)
        log(f"Elapsed: {elapsed} ns")

        print(f"{sample_size},{elapsed}", file=f)
        sample_size *= 2

    f.close() 

