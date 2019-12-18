#! /usr/bin/env python
import argparse
import sys
import csv
import pprint
import os

def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_taxonomy_tsv')
    p.add_argument('filenames')
    p.add_argument('accs')
    args = p.parse_args()

    fp = open(args.gtdb_taxonomy_tsv)
    gtdb_r = csv.DictReader(fp, delimiter='\t', fieldnames=['ident', 'gtdb_tax'])
    gtdb_ident_to_tax = {}
    gtdb_ident_to_filename = {}
    for row in gtdb_r:
        ident = row['ident']
        assert ident not in gtdb_ident_to_tax # no dups?
        gtdb_ident_to_tax[ident] = row['gtdb_tax']

        if ident.startswith('GB_') or ident.startswith('RS_'):
            filename = 'database/' + ident[3:] + '_genomic.fna.gz'
        elif ident.startswith('UBA'):
            filename = 'database/' + ident + '_genomic.fna.gz'
        else:
            assert 0, ident
        gtdb_ident_to_filename[ident] = filename

    print('loaded {} tax entries'.format(len(gtdb_ident_to_tax)))

    fp = open(args.filenames)
    fp2 = open(args.accs)
    filenames_to_acc = {}
    for (filename, acc) in zip(fp, fp2):
        filename = filename.strip()
        acc = acc.strip()
        filenames_to_acc[filename] = acc

    print('loaded {} filenames & accs'.format(len(filenames_to_acc)))

    gtdb_ident_to_acc = {}
    for ident, filename in gtdb_ident_to_filename.items():
        acc = filenames_to_acc[filename]
        gtdb_ident_to_acc[ident] = acc

    # output new spreadsheet
    fp = open('gtdb_files/gtdb-new-tax.tsv', 'wt')
    w = csv.writer(fp, delimiter='\t')
    w.writerow(['identifier', 'taxonomy', 'filename', 'first_accession'])
    for ident in gtdb_ident_to_tax:
        tax = gtdb_ident_to_tax[ident]
        filename = gtdb_ident_to_filename[ident]
        acc = gtdb_ident_to_acc[ident]

        w.writerow([ident, tax, filename, acc])

    print('created gtdb_files/gtdb-new-tax.tsv')

    fp = open('gtdb_files/gtdb-lineages.csv', 'wt')
    w = csv.writer(fp)
    w.writerow('accession,gtdb_id,superkingdom,phylum,class,order,family,genus,species'.split(','))
    for ident, tax in gtdb_ident_to_tax.items():
        acc = gtdb_ident_to_acc[ident]
        acc = acc.split('.')[0]
        row = [acc, ident] + tax.split(';')

        w.writerow(row)

    print('created gtdb-lineages.csv')
    
if __name__ == '__main__':
    main()
