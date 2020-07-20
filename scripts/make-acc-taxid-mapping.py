#! /usr/bin/env python
"""
Take an NCBI 'assembly_summary.txt' file and print to
stdout the accessions and their associated taxids.
"""

import argparse
import csv
import gzip


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--strip-version', action='store_true')
    p.add_argument('summary')
    args = p.parse_args()

    with gzip.open(args.summary, 'rt') as fp:
        fp.readline() # skip first line
        fp.read(2) # skip initial comment in header
        data = csv.DictReader(fp, delimiter='\t')
        for row in data:
            accession = row['assembly_accession']
            if args.strip_version:
                # --split-identifiers in `sourmash lca index` doesn't behave
                # well with version, so remove it
                accession = accession.split('.')[0]
            taxid = row['taxid']
            print(f'{accession},{taxid}')


if __name__ == '__main__':
    main()
