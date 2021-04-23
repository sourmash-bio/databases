import os
import sys
import argparse

#import sourmash
import pandas as pd


def process_infofile(inF):
    # Read and modify single taxonomy csv
    info_csv = pd.read_csv(inF, sep = "\t", names=['ident', 'gtdb_tax'])
    # remove RS_, GB_ that may exist in front of accession
    info_csv["ident"] = info_csv["ident"].str.replace("RS_", "").str.replace("GB_", "") #.str.rsplit(".", 1, expand=True)[0]
    #info_csv["accession"] = info_csv["ident"].str.rsplit(".", 1, expand=True)[0]
    info_csv[["superkingdom","phylum","class","order","family","genus","species"]] = info_csv["gtdb_tax"].str.split(pat=";", expand=True)
    info_csv.drop(columns=["gtdb_tax"], inplace=True)
    print(f"loaded {len(info_csv)} tax entries from {inF}")
    return info_csv



#for ident, tax in gtdb_ident_to_tax.items():
#
#        acc = ident.split('.')[0]
#        filename = ident + '_genomic.fna.gz'
#        row = [acc,filename] + tax.split(';')
#
#        w.writerow(row)




def main(args):
    # process each taxonomy file, then concat to a single dataframe
    frames = [process_infofile(f) for f in args.taxonomy_files]
    full_info = pd.concat(frames)
    full_info.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--taxonomy-files", nargs="+", help="gtdb taxonomy files to merge (bac120_taxonomy, ar122_taxonomy)")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
