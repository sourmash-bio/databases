#! /usr/bin/env python
"""
Retrieve all the accessions from a list of genomes.
"""
import sys
import argparse
import gzip


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genome_listfile', nargs='+')
    args = p.parse_args()

    for filename in args.genome_listfile:
        with open(filename, 'rt') as fp:
            for genomefile in fp:
                line = next(iter(gzip.open(genomefile.strip())))
                line = line.split()[0]
                line = line[1:]
                
                print(line.decode('utf-8'))
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
