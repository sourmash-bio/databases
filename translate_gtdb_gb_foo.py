#! /usr/bin/env python
import argparse
import sys
import csv
import pprint
import os

DOMAINS = ['bacteria', 'fungi', 'archaea', 'viral']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_taxonomy_tsv')
    p.add_argument('gb_assembly_summary')
    args = p.parse_args()

    ###

    filenames_to_acc = {}
    for domain in DOMAINS:
        with open('domain-{}.txt'.format(domain), 'rt') as fp:
            d_filenames = [ os.path.basename(x.strip()) for x in fp ]
        with open('domain-{}.acc.txt'.format(domain), 'rt') as fp:
            d_accs = [ x.strip() for x in fp ]

        d = zip(d_filenames, d_accs)
        e = dict(d)
        filenames_to_acc.update(e)

    print('loaded filenames to acc:', len(filenames_to_acc))

    ###

    fp = open(args.gb_assembly_summary, 'rt')
    fieldnames = "# assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material"[2:].split('\t')
    gb_r = csv.DictReader(fp, delimiter='\t', fieldnames=fieldnames)

    next(gb_r)
    next(gb_r)

    print('loading genbank assembly summary...')
    gb_asmacc_to_rs = {}
    gb_rs_to_asmacc = {}
    asmacc_to_filename = {}
    for row in gb_r:
        asmacc = row['assembly_accession']
        rs = row['gbrs_paired_asm']
        gb_asmacc_to_rs[asmacc] = rs

        filename = os.path.basename(row['ftp_path']) + '_genomic.fna.gz'
        asmacc_to_filename[asmacc] = filename

        if rs.strip():
            gb_rs_to_asmacc[rs] = asmacc
        
    print('loaded GB:', len(gb_asmacc_to_rs))
    print('entries with RefSeq mapping:', len(gb_rs_to_asmacc))
    
    fp = open(args.gtdb_taxonomy_tsv)
    gtdb_r = csv.DictReader(fp, delimiter='\t', fieldnames=['ident', 'gtdb_tax'])
    gtdb_ident_to_tax = {}
    for row in gtdb_r:
        assert row['ident'] not in gtdb_ident_to_tax # no dups?
        gtdb_ident_to_tax[row['ident']] = row['gtdb_tax']

    print('loaded GTDB:', len(gtdb_ident_to_tax))

    found_rs = 0
    missing_rs = 0
    non_gb = 0
    gtdb_ident_to_asmacc = {}
    for ident in gtdb_ident_to_tax:
        if ident.startswith('RS_'):
            rs_id = ident[3:]
            if rs_id in gb_rs_to_asmacc:
                found_rs += 1
                gtdb_ident_to_asmacc[ident] = gb_rs_to_asmacc[rs_id]
            else:
                missing_rs += 1
            
        elif ident.startswith('GB_'):
            gtdb_ident_to_asmacc[ident] = ident[3:]
        else:
            non_gb += 1

    asmacc_to_sigacc = {}
    for asmacc, filename in asmacc_to_filename.items():
        sigacc = filenames_to_acc.get(filename)
        if sigacc:
            asmacc_to_sigacc[asmacc] = sigacc

    gtdb_ident_to_sigacc = {}
    x = 0
    y = 0
    for ident, asmacc in gtdb_ident_to_asmacc.items():
        sigacc = asmacc_to_sigacc.get(asmacc)
        if sigacc:
            x += 1
            gtdb_ident_to_sigacc[ident] = sigacc
        else: y += 1

    print(x, y)

    print('')
    print('total in GTDB tax:', len(gtdb_ident_to_tax))
    print('found RefSeq -> GenBank mapping:', found_rs)
    print('missing Refseq -> GenBank mapping:', missing_rs)
    print('non GenBank/RefSeq accessions in GTDB:', non_gb)

    x = 0
    y = 0
    for ident in gtdb_ident_to_tax:
        if ident in gtdb_ident_to_sigacc:
            x += 1
        else:
            y += 1
    print('summary:\ntotal GTDB tax assignments:', x)
    print('num not matched to sourmash sig accession id:', y)

    fp = open('gtdb-lineages.csv', 'wt')
    w = csv.writer(fp)
    w.writerow('accession,gtdb_id,superkingdom,phylum,class,order,family,genus,species'.split(','))

    for ident, tax in gtdb_ident_to_tax.items():
        sigacc = gtdb_ident_to_sigacc.get(ident)
        if sigacc:
            sigacc = sigacc.split('.')[0]
            tax = tax.split(';')
            row = [sigacc, ident] + tax

            w.writerow(row)


if __name__ == '__main__':
    sys.exit(main())
