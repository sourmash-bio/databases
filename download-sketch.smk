'''
Snakefile for downloading and sketching genomes needed for building sourmash databases
'''

import csv
import pandas as pd
import numpy as np

out_dir = config.get('output_dir', "output.databases")
logs_dir = os.path.join(out_dir, 'logs')

metadata_file = config['metadata'] 
metadata = pd.read_csv(metadata)
basename = config['basename']

ACCESSIONS = metadata['ident']


# checkpoint function to read the fromfile csv and handle fasta downloads if necessary
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern
    
    def get_filenames(self, basename=None, moltype=None):
        df = pd.read_csv(f"{out_dir}/{basename}.fromfile.csv")
        filename_col = 'genome_filename'
        if moltype == "protein":
            filename_col = 'protein_filename'
        # filter df to get non-empties in relevant column
        fastas = df[filename_col][df[filename_col].notnull()].tolist()
        return fastas

    def __call__(self, w):
        global checkpoints
        # wait for the results of 'check_fromfile'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_fromfile.get(**w)
        #expand the pattern
        fastas = self.get_filenames(**w)

        pattern = expand(self.pattern, fn =fastas,**w)
        return pattern

rule all:
    input:
        expand(os.path.join(out_dir, f"{basename}.{{moltype}}.zip"), moltype = ['dna', 'protein']),

### Rules for Assemblies:
# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = 'genbank/info/{acc}.info.csv'
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=90,
        time=60,
        partition="low2",
    conda: "conf/env/ncbi.yml"
    shell:
        """
        python -Werror genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
        """

# download actual genomes!
rule download_matching_genome_wc:
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        genome = protected("genbank/genomes/{acc}_genomic.fna.gz")
    threads: 1
    resources:
        mem_mb=3000,
        disk_mb=5000,
        runtime=60,
        time=90,
        partition="low2",
    # conda: "conf/env/genomes.yml"
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...",
                    file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                           file=sys.stderr)


# download proteome (when it exists)
rule download_matching_proteome_wc:
    input:
        csvfile = ancient('genbank/info/{acc}.info.csv')
    output:
        proteome = protected("genbank/proteomes/{acc}_protein.faa.gz")
    threads: 1
    resources:
         mem_mb=3000,
         runtime=60,
         time=90,
         partition="low2", # bml
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['protein_url']
            name = row['ncbi_tax_name']

            print(f"downloading proteome for acc {acc}/{name} from NCBI...",
                file=sys.stderr)
            with open(output.proteome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.proteome}",
                          file=sys.stderr)


localrules: build_fromfile_from_assemblyinfo
rule build_fromfile_from_assemblyinfo:
    input: 
        info=expand("genbank/info/{acc}.info.csv", acc=ACCESSIONS)
    output:
        csv = os.path.join(out_dir, "{basename}.fromfile.csv")
    threads: 1
    run:
        with open(str(output.csv), "w") as outF:
            header = ["name","genome_filename","protein_filename"]
            outF.write(','.join(header) + '\n')
            for inp in input:
                with open(str(inp)) as csvfile:
                    csv_reader = csv.DictReader(csvfile)# .# has fields:  ["acc", "genome_url", "protein_url", "assembly_report_url", "ncbi_tax_name", "ncbi_taxid"]
                    rows = list(csv_reader)
                    assert len(rows) == 1
                    row = rows[0]
                    acc = row['acc']
                    name = acc + " " + row['ncbi_tax_name']
                    gf = f"genbank/genomes/{acc}_genomic.fna.gz"
                    pf= ""
                    if row["protein_url"]:
                        # do we want to add prodigal to go nucl --> protein where we don't yet have protein fastas to download?
                        pf = f"genbank/proteomes/{acc}_protein.faa.gz"
                    outF.write(f"{name},{gf},{pf}\n")

    
 # Define the checkpoint function that allows us to read the fromfile.csv
checkpoint check_fromfile:
    input: os.path.join(out_dir, f"{basename}.fromfile.csv"),
    output: touch(os.path.join(out_dir,".check_fromfile"))


# could enable more here
paramD = {"dna": "dna,k=21,k=31,51,scaled=1000,abund", "protein": "protein,k=7,k=10,scaled=200,abund"}
rule sketch_fromfile:
    input: 
        fromfile=os.path.join(out_dir, "{basename}.fromfile.csv"),
        fastas=ancient(Checkpoint_MakePattern("{fn}")),
    output: os.path.join(out_dir, "{basename}.{moltype}.zip")
    params:
        lambda w: paramD[w.moltype]
    threads: 1
    resources:
        mem_mb=6000,
        runtime=4000,
        time=4000,
        partition="high2",
    conda: "conf/env/sourmash.yml"
    log:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.log")
    benchmark:  os.path.join(logs_dir, "sketch", "{basename}.{moltype}.benchmark")
    shell:
        """
        sourmash sketch fromfile {input.fromfile} -p {params} -o {output} --report-duplicated --ignore-missing 2> {log}
        """

