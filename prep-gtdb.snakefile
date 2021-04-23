"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  prep-gtdb.snakefile -n
"""

import os
import sys
import pandas as pd

configfile: "gtdb-r202.yml"

out_dir = config['output_dir']
#sig_dir = config['sig_dir']
logs_dir = os.path.join(out_dir, "logs")
data_dir = os.path.join(out_dir, "inputs")

#taxon_info = config["taxa"]
# basename-r{release}
basename = config["basename"] + "-r" + str(config["release"])

# check params are in the right format
for alpha, info in config["alphabet_info"].items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        config["alphabet_info"][alpha]["scaled"] = [scaled]
    if not isinstance(ksize, list):
        config["alphabet_info"][alpha]["ksize"] = [ksize]

# ctb checkpoint code to specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def check_for_sigs(self):
        global tax_info
        tax_info = pd.read_csv(f'{out_dir}/{basename}.taxonomy.csv', header=0)
        accessions = tax_info['ident']
        #accessions = tax_info['accession']
        if config["build_genomic"]:
            tax_info["genomic_filename"] = tax_info["ident"] + config["genomic"]["fasta_suffix"]
            # to do:
              #- check for files that exist; make dictionary of *just those that don't exist yet*
              #- use that dictionary to download genomes + sketch --> wort dir OR new dir?
              #


        if config["build_protein"]:
            tax_info["protein_filename"] = tax_info["ident"] + config["protein"]["fasta_suffix"]
        tax_info.set_index("ident", inplace=True)
        return accessions

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_csv.get(**w)

        accs = self.check_for_sigs()

        pattern = expand(self.pattern, accession=accs, **w)
        return pattern

rule all:
    #input: os.path.join(out_dir, f"{basename}.siglist.txt")
    #input: expand(os.path.join(data_dir, "{taxfile}"), taxfile=config["taxonomy_files"].keys())
    #input: os.path.join(out_dir, f"{basename}.taxonomy.csv") 
    input: os.path.join(out_dir, f"{basename}.genomic.siglist.txt")

localrules: download_csvs
rule download_csvs:
    message: "Download GTDB CSVs"
    output:
        expand(os.path.join(data_dir, "{taxfile}"), taxfile=config["taxonomy_files"].keys())
    params:
        taxinfo = config["taxonomy_files"],
    run:
        for outfile, link in params.taxinfo.items():
            full_outfile = os.path.join(data_dir, outfile)
            shell("wget -O {full_outfile} {link}")

localrules: make_taxonomy_csv
rule make_taxonomy_csv:
    message: "Merge and Format CSVs into a single taxonomy csv"
    input: rules.download_csvs.output
    output: os.path.join(out_dir, f"{basename}.taxonomy.csv")
    log: os.path.join(logs_dir, "make_taxonomy_csv", f"{basename}.make_taxonomy_csv.log")
    threads: 1
    shell:
        """
        python make-gtdb-taxonomy.py --taxonomy-files {input} --output-csv {output} > {log} 2>&1
        """

# make the checkpoint work
localrules: check_csv
checkpoint check_csv:
    input: os.path.join(out_dir, f"{basename}.taxonomy.csv")
    output: touch(f"{out_dir}/.check_csv")


'''
# download genomes and build signatures if they don't already exist?
rule ncbi_datasets_download:
    input: 
        check_csv=f"{out_dir}/.check_csv",
        tool="scripts/ncbi-datasets"
    output:
        genomic=protected(os.path.join(out_dir, "data/{accession}.fna.gz")),
        #protein=protected(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"))
    params:
        tmp= lambda w: os.path.join(out_dir, w.accession) + ".zip"
    threads: 1
    resources:
        mem_mb= 2000,
        runtime= 60
    log: os.path.join(logs_dir, "ncbi_datasets", "{accession}.download")
    shell:
        """
        echo "scripts/ncbi-datasets download genome accession {wildcards.accession} --exclude-rna --exclude-gff3 -f {params.tmp}" > {log} 
        scripts/ncbi-datasets download genome accession {wildcards.accession} --exclude-rna --exclude-gff3 -f {params.tmp}
        unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*.fna | gzip -9 > {output.genomic}
        unzip {params.tmp} -d nd_{wildcards.accession}
        echo "fna filenames: " >> {log}
        ls -1 nd_{wildcards.accession}/ncbi_dataset/data/{wildcards.accession}/*.fna >> {log}
        rm -rf nd_{wildcards.accession} 
        rm -rf {params.tmp} 
        """
        #unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*.faa | gzip -9 > {output.protein}

def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule sourmash_sketch_dna:
    input:
        ancient(os.path.join(data_dir, "{accession}.fna.gz")),
    output:
        os.path.join(config["genomic"]["sig_dir"], "{accession}.sig"),
    params:
        sketch_params=make_param_str(config["alphabet_info"]["nucleotide"]["ksize"], config["alphabet_info"]["nucleotide"]["scaled"]),
        signame = lambda w: tax_info.at[w.accession, "signame"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """
'''


localrules: signames_to_file_genomic
rule signames_to_file_genomic:
    input: 
        csv=f"{out_dir}/.check_csv",
        sigs=Checkpoint_MakePattern(os.path.join(config["genomic"]["sig_dir"], "{accession}.sig")),
    output: os.path.join(out_dir, f"{basename}.genomic.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                #outF.write(str(inF) + "\n")
                outF.write(full_filename + "\n")

