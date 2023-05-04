import os
import sys
import argparse
import pandas as pd


def process_infofile(inF):
    # Read and modify single taxonomy csv
    info_csv = pd.read_csv(inF, sep = "\t", names=['ident', 'gtdb_tax'])
    # remove RS_, GB_ that may exist in front of accession. keep full ident so we can download properly if sigfile doesn't already exist
    info_csv["ident"] = info_csv["ident"].str.replace("RS_", "").str.replace("GB_", "") #.str.rsplit(".", 1, expand=True)[0]
    #info_csv["accession"] = info_csv["ident"].str.rsplit(".", 1, expand=True)[0]
    info_csv[["superkingdom","phylum","class","order","family","genus","species"]] = info_csv["gtdb_tax"].str.split(pat=";", expand=True)
    info_csv.drop(columns=["gtdb_tax"], inplace=True)
    print(f"loaded {len(info_csv)} tax entries from {inF}")
    return info_csv
        
def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--metadata-files", nargs="+", help="gtdb metadata files")
    p.add_argument("--output-prefix",  help="output csv files prefix")
    args = p.parse_args()
    return main(args)


def main(args):
    
    current_release = f"{args.output_prefix}.current.csv"
    diff_release = f"{args.output_prefix}.diff.csv"
    
    metadata_dfs = []
    for metadata_file in args.metadata_files:
        metadata_info = pd.read_csv(metadata_file, header=0, low_memory=False, sep = "\t") 
        filtered_metadata_info = metadata_info[["accession", "gtdb_representative", "gtdb_taxonomy"]]
        filtered_metadata_info["accession"] = metadata_info["accession"].str.replace("RS_", "").str.replace("GB_", "") #.str.rsplit(".", 1, expand=True)[0]
        metadata_dfs.append(filtered_metadata_info)
        

    # Writing thew new metadata csv file
    metadata_info = pd.concat(metadata_dfs)
    metadata_info[["superkingdom","phylum","class","order","family","genus","species"]] = metadata_info["gtdb_taxonomy"].str.split(pat=";", expand=True)
    metadata_info.drop(columns=["gtdb_taxonomy"], inplace=True)
    metadata_info.rename(columns={"accession":"ident"}, inplace=True)
    metadata_info.to_csv(current_release, sep = '\t', index=False)
    
    
    # Load the old version's metadata csv file
    old_metadata_info = pd.read_csv("/group/ctbrowngrp/sourmash-db/gtdb-rs207/gtdb-rs207.taxonomy.csv", sep = ",")
    
    # Create a new df for difference in accesisons between two csv files
    new_df = metadata_info[~metadata_info.ident.isin(old_metadata_info.ident)]
    new_df.to_csv(diff_release, sep = ',', index=False)


if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
