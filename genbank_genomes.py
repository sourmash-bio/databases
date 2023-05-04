#! /usr/bin/env python
"""
Modified from:
 - https://raw.githubusercontent.com/bluegenes/2021-rank-compare/main/genbank_genomes.py

Compared with genome-grist version, this version:
  - also includes protein faa download link
  - checks that each url exists before using it for download
  - includes NBCI tax_name and taxid
"""
import sys
import argparse
import urllib.request
import csv

from lxml import etree

def url_exists(url):
    """
    Given a URL, check if the URL returns a 200 status code (i.e. exists) when accessed via HTTP.
    """
    try:
        response = urllib.request.urlopen(url)
        return response.getcode() == 200
    except urllib.error.HTTPError as e:
        return False

def url_for_accession(accession):
    db, acc = accession.strip().split("_")
    if '.' in acc:
        number, version = acc.split(".")
    else:
        number, version = acc, '1'
    number = "/".join([number[p : p + 3] for p in range(0, len(number), 3)])
    url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"
    print(f"opening directory: {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as response:
        all_names = response.read()

    print("done!", file=sys.stderr)

    all_names = all_names.decode("utf-8")

    full_name = None
    for line in all_names.splitlines():
        if line.startswith(f'<a href='):
            name=line.split('"')[1][:-1]
            db_, acc_, *_ = name.split("_")
            if db_ == db and acc_.startswith(acc):
                full_name = name
                break

    if full_name is None:
        return None
    else:
        url = "htt" + url[3:]
        urls = [f"{url}/{full_name}/{full_name}_genomic.fna.gz",
                f"{url}/{full_name}/{full_name}_protein.faa.gz",
                f"{url}/{full_name}/{full_name}_assembly_report.txt"]
        existing_urls = []
        for url in urls:
            if url_exists(url):
                existing_urls.append(url)
            else:
                existing_urls.append("")
        return(tuple(existing_urls))

def get_taxid_from_assembly_report(url):
    print(f"opening assembly report: {url}", file=sys.stderr)
    with urllib.request.urlopen(url) as response:
        content = response.read()
    print("done!", file=sys.stderr)

    content = content.decode("utf-8").splitlines()
    for line in content:
        if "Taxid:" in line:
            line = line.strip()
            pos = line.find("Taxid:")
            assert pos >= 0
            pos += len("Taxid:")
            taxid = line[pos:]
            taxid = taxid.strip()
            return taxid

    assert 0


def get_tax_name_for_taxid(taxid):
    tax_url = (
        f"https://www.ncbi.nlm.nih.gov/taxonomy/?term={taxid}&report=taxon&format=text"
    )
    print(f"opening tax url: {tax_url}", file=sys.stderr)
    with urllib.request.urlopen(tax_url) as response:
        content = response.read()

    print("done!", file=sys.stderr)

    root = etree.fromstring(content)
    notags = etree.tostring(root).decode("utf-8")
    if notags.startswith("<pre>"):
        notags = notags[5:]
    if notags.endswith("</pre>"):
        notags = notags[:-6]
    notags = notags.strip()

    return notags


def main(args):
    fieldnames = ["acc", "genome_url", "protein_url", "assembly_report_url", "ncbi_tax_name", "ncbi_taxid"]
    fp=None
    if args.output:
        fp = open(args.output, "wt")
        w = csv.DictWriter(fp, fieldnames=fieldnames)
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    w.writeheader()

    acc = args.accession

    genome_url, protein_url, assembly_report_url = url_for_accession(acc)
    taxid = get_taxid_from_assembly_report(assembly_report_url)
    tax_name = get_tax_name_for_taxid(taxid)

    d = dict(
        acc=acc,
        genome_url=genome_url,
        protein_url=protein_url,
        assembly_report_url=assembly_report_url,
        ncbi_tax_name=tax_name,
        ncbi_taxid=taxid,
    )

    w.writerow(d)
    print(f"retrieved for {acc} - {tax_name}", file=sys.stderr)
    if fp:
        fp.close()

    return 0

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("accession")
    p.add_argument("-o", "--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
