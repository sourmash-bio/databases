#! /usr/bin/env python
import argparse
import sys
import csv
import pprint
import os

def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_taxonomy_tsv')
    args = p.parse_args()

    fp = open(args.gtdb_taxonomy_tsv)
    gtdb_r = csv.DictReader(fp, delimiter='\t', fieldnames=['ident', 'gtdb_tax'])
    gtdb_ident_to_tax = {}
    for row in gtdb_r:
        ident = row['ident']
        assert ident not in gtdb_ident_to_tax # no dups?

        if ident.startswith('GB_') or ident.startswith('RS_'):
            ident = ident[3:]
        elif ident.startswith('UBA'):
            pass
        else:
            assert 0, ident

        gtdb_ident_to_tax[ident] = row['gtdb_tax']

    print('loaded {} tax entries'.format(len(gtdb_ident_to_tax)))

    fp = open('gtdb_files/gtdb-lineages.csv', 'wt')
    w = csv.writer(fp)
    w.writerow('accession,gtdb_id,superkingdom,phylum,class,order,family,genus,species'.split(','))
    for ident, tax in gtdb_ident_to_tax.items():
        row = [ident] + tax.split(';')

        w.writerow(row)

    print('created gtdb_files/gtdb-lineages.csv')
    
if __name__ == '__main__':
    main()
