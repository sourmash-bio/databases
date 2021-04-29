"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  prep-gtdb.snakefile --group-components sketch=5000 -n
"""
#- check for files that exist; make dictionary of *just those that don't exist yet*
#- use that dictionary to download genomes + sketch --> local dir?
#- write siglist with correct path for each sigfile

import os
import sys
import pandas as pd

configfile: "gtdb-r202.yml"

out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
data_dir = os.path.join(out_dir, "inputs")

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

    

nucl_siglists, prot_siglists = [],[]
# genomic zipfiles
if config.get("build_genomic", False):
    nucl_siglists.append("genomic")
    if config.get("build_representative_sets", False):
        nucl_siglists.append("genomic-reps")
# protein zipfiles
if config.get("build_protein", False):
    prot_siglists.append("protein")
    if config.get("build_representative_sets", False):
        prot_siglists.append("protein-reps")

siglist_types = nucl_siglists + prot_siglists
# temp hack for protein rep sets
siglist_types = ["protein-reps"]

# set signature directories
existing_genomic_sigdir = config["genomic"]["sig_dir"]
new_genomic_sigdir = os.path.join(out_dir, "genomic", "signatures")

existing_protein_sigdir = config["protein"]["sig_dir"]
new_protein_sigdir = os.path.join(out_dir, "protein", "signatures")
prodigal_sigdir = os.path.join(data_dir, "prodigal", "signatures")


def sigfile_exists(sigpath):
    # does this sigfile exist already? 
    return os.path.isfile(sigpath)

def update_sigpath(sigfile, sigdir="", new_sigdir=""):
    sigpath = os.path.join(sigdir, sigfile)
    if os.path.isfile(sigpath):
        return sigpath
    else:
        new_sigpath = os.path.join(new_sigdir, sigfile)
        return new_sigpath


class Checkpoint_MakePattern:

    def __init__(self, pattern):
        self.pattern = pattern

    def build_taxinfo(self, basename=None, input_type=None):
        tax_info = pd.read_csv(f"{out_dir}/{basename}.taxonomy.csv", header=0)
        #For metadata info, pandas needs low_memory=False(DtypeWarning: Columns (61,65,74,82,83) have mixed types.Specify dtype option on import or set low_memory=False)
        metadata_info = pd.read_csv(f"{out_dir}/{basename}.metadata.csv.gz", header=0, low_memory=False) 
        representative_accessions = metadata_info[metadata_info["gtdb_representative"] == "t"]["accession"].str.replace("GB_", "").str.replace("RS_", "")
        tax_info["signame"] = tax_info["ident"] + " " + tax_info["species"].str.replace("s__", "")
        # build sigfile basenames
        tax_info["genomic_sigfile"] = tax_info["ident"] + ".sig"
        tax_info["protein_sigfile"] = tax_info["ident"] + ".sig"
        # set accessions as DataFrame Index
        tax_info.set_index("ident", inplace=True)
        return tax_info, representative_accessions

    def find_sigs(self, basename=None, input_type=None):
        # Update sigfile paths for accessions we need to download and sketch (genomic, protein)
        
        if "genomic" in input_type:
            print("finding genomic sigs...")
            sigpath_name = "genomic_sigfile"
            existing_sigdir = existing_genomic_sigdir
            new_sigdir = new_genomic_sigdir
        elif "protein" in input_type:
            print("finding protein sigs...")
            sigpath_name = "protein_sigfile"
            existing_sigdir = existing_protein_sigdir
            new_sigdir = new_protein_sigdir
        
        # first, update sigpath -- is it in the wort dir, or our local dir?
        tax_info[sigpath_name] = tax_info[sigpath_name].apply(update_sigpath, sigdir=existing_sigdir, new_sigdir=new_sigdir)
        # just for nice output, figure out the number of sigfiles that still need to be sketched + print 
        sketchD = tax_info.loc[~tax_info["genomic_sigfile"].apply(sigfile_exists)]
        print(f"found {len(sketchD)}accessions that need to be downloaded and sketched for {input_type} databases")
        # return the path for all sigfiles needed for this database
        if "reps" in input_type:
            representatives = tax_info.loc[tax_info.index.isin(representative_accessions)]
            return representatives[sigpath_name]
        else:
            return tax_info[sigpath_name]

    
    def __call__(self, w):
        global checkpoints
        global tax_info
        global representative_accessions
        tax_info = None
        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_csv.get(**w)

        if not tax_info:
            tax_info, representative_accessions = self.build_taxinfo(**w)
        
        sigfiles = self.find_sigs(**w)
        pattern = expand(self.pattern, sigfile=sigfiles, **w)
        
        return pattern


class Check_Protein_Sigs:
    def __init__(self, pattern):
        self.pattern = pattern

    def update_sigpaths(self, basename=basename, input_type=None):
        with open(f"{out_dir}/{basename}.{input_type}.prodigal-siglist.txt", "rt") as fp:
            accessions = [ os.path.basename(x).rstrip().rsplit(".sig")[0] for x in fp ]
            for accession in accessions:
                tax_info.at[accession, "protein_sigfile"] = f"{prodigal_sigdir}/{accession}.prodigal.sig"
            if input_type == "protein":
                return tax_info["protein_sigfile"]
            elif input_type == "protein-reps":
                # find representative protein files
                representative_protein = tax_info.loc[tax_info.index.isin(representative_accessions)]
                return representative_protein["protein_sigfile"]
            else:
                print(f"Unknown input type {input_type}!")
                sys.exit(-1)
            
    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_proteins'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_proteins.get(**w)
        prot_sigfiles = self.update_sigpaths(**w)

        pattern = expand(self.pattern, sigfile=prot_sigfiles, **w)
        return pattern



rule all:
    input: 
        expand(os.path.join(out_dir, "databases", "{basename}.{input_type}.zip"), basename=basename, input_type=siglist_types),
        expand(os.path.join(out_dir, "ksize_databases", "{basename}.{input_type}.k{ksize}.zip"), basename=basename, input_type=nucl_siglists, ksize=config["alphabet_info"]["nucleotide"]["ksize"])

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

localrules: download_metadata_and_combine
rule download_metadata_and_combine:
    message: "Download GTDB Metadata"
    output: os.path.join(out_dir, f"{basename}.metadata.csv.gz")
    params:
        outdir = data_dir + "/", # wget --> tar redirect fails without trailing "/"
        metadata_info = config["metadata_files"],
    run:
        for link in params.metadata_info.values():
            shell("wget {link} -O - | tar  -xzvf - -C {params.outdir}")
        frames = []
        for f in params.metadata_info.keys():
            tsv_file = params.outdir + f
            frames+=[pd.read_csv(tsv_file, sep="\t", header=0, low_memory=False)]
        full_info = pd.concat(frames)
        full_info.to_csv(str(output), index=False)
 
# make the checkpoint work
localrules: check_csv
checkpoint check_csv:
    input: 
        taxonomy=os.path.join(out_dir, "{basename}.taxonomy.csv"),
        metadata=os.path.join(out_dir, "{basename}.metadata.csv.gz")
    output: 
        touch(os.path.join(out_dir, ".{basename}.{input_type}.check_csv"))

# from https://github.com/dib-lab/sourmash_databases/blob/6e93e871f2e955853b23e54c12b4fc42fb26ef1c/Snakefile.assembly
def url_for_accession(accession, protein=False):
    fasta_ext = "genomic.fna.gz"
    if protein:
        fasta_ext = "protein.faa.gz"
    db, acc = accession.split("_")
    number, version = acc.split(".")
    number = "/".join([number[pos:pos + 3] for pos in range(0, len(number), 3)])
    url = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{number}"

    from subprocess import CalledProcessError
    try:
        all_names = shell(f"curl -s -l {url}/", read=True).split('\n')
    except CalledProcessError as e:
        #raise ValueError(f"Can't find URL for {accession}, tried {url}")
        return "" #ignore and handle later!

    full_name = None
    for name in all_names:
        db_, acc_, *_ = name.split("_")
        if db_ == db and acc == acc_:
            full_name = name
            break

    url = "https" + url[3:]
    return f"{url}/{full_name}/{full_name}_{fasta_ext}"


def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule stream_and_sketch_genomic:
    output:
        os.path.join(new_genomic_sigdir, "{accession}.sig"),
    params:
        sketch_params=make_param_str(config["alphabet_info"]["nucleotide"]["ksize"], config["alphabet_info"]["nucleotide"]["scaled"]),
        signame = lambda w: tax_info.at[w.accession, "signame"],
        url_path = lambda w: url_for_accession(w.accession),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_genomic", "{accession}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_genomic", "{accession}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} <(curl -s {params.url_path} | zcat) 2> {log}
        """

def make_protein_param_str(alphabets=["protein", "dayhoff", "hp"]):
    param_str = ""
    for alpha in alphabets:
        info = config["alphabet_info"][alpha]
        ksizes = info["ksize"]
        scaled = info["scaled"]
        ks = [ f'k={k}' for k in ksizes ]
        ks = ",".join(ks)
        scaled = min(scaled) #take minimum value of scaled list
        param_str += f" -p {alpha},{ks},scaled={scaled},abund "
    return param_str

rule stream_and_sketch_protein:
    output:
        os.path.join(new_protein_sigdir, "{accession}.sig"),
    params:
        sketch_params = make_protein_param_str(),
        signame = lambda w: tax_info.at[w.accession, "signame"],
        url_path = lambda w: url_for_accession(w.accession, protein=True),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_protein", "{accession}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_protein", "{accession}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch protein {params.sketch_params} -o {output} --name {params.signame:q} <(curl -s {params.url_path} | zcat) 2> {log}
        touch {output} 
        """

rule check_protein_sigfiles_for_empties:
    input:
        csv=os.path.join(out_dir, ".{basename}.{input_type}.check_csv"),
        sigs=ancient(Checkpoint_MakePattern("{sigfile}")),
    output: os.path.join(out_dir, "{basename}.{input_type}.prodigal-siglist.txt")
    wildcard_constraints:
        input_type="protein|protein-reps"
    run:
        with open(str(output), 'w') as out:
            for sig in input.sigs:
                sig_str = str(sig)
                if os.path.getsize(sig_str) == 0:
                    out.write(f"{sig_str}\n")


# make the checkpoint work
localrules: check_proteins
checkpoint check_proteins:
    input:
        os.path.join(out_dir, "{basename}.{input_type}.prodigal-siglist.txt")
    wildcard_constraints:
        input_type="protein|protein-reps"
    output: touch(os.path.join(out_dir, ".{basename}.{input_type}.check_proteins"))
          

rule download_genomes_for_failed_protein_sigs:
    output: temp(os.path.join(data_dir, "{accession}_genomic.fna")) # output file marked as temp is deleted after all rules that use it as an input are completed
    params:
        url_path = lambda w: url_for_accession(w.accession),
    log: os.path.join(logs_dir, "download_genomes_for_failed_protein_sigs", "{accession}.download")
    shell:
        """
        curl -s {params.url_path} | gunzip > {output} 2> {log}
        """

rule prodigal_genomes_for_failed_protein_sigs:
    input: os.path.join(data_dir, "{accession}_genomic.fna")
    output: 
        proteins=temp(os.path.join(data_dir, "{accession}.prodigal.faa")),
        genes=temp(os.path.join(data_dir, "{accession}.prodigal.fna"))
    conda: "envs/prodigal-env.yml"
    log: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.log")
    benchmark: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.benchmark")
    shell:
        """
         prodigal -i {input} -a {output.proteins} -o {output.genes} 2> {log}
        """

rule sketch_prodigal_genomes_for_failed_protein_downloads:
    input: os.path.join(data_dir, "{accession}.prodigal.faa")
    output: os.path.join(prodigal_sigdir,  "{accession}.prodigal.sig"),
    params:
        sketch_params = make_protein_param_str(),
        signame = lambda w: tax_info.at[w.accession, "signame"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    group: "sketch"
    log: os.path.join(logs_dir, "sourmash_sketch_protein", "{accession}.prodigal-sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_protein", "{accession}.prodigal-sketch.benchmark")
    conda: "envs/sourmash4.yml"
    shell:
        """
        sourmash sketch protein {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """

localrules: genomic_signames_to_file
rule genomic_signames_to_file:
    input: 
        csv=os.path.join(out_dir, ".{basename}.{input_type}.check_csv"),
        sigs=ancient(Checkpoint_MakePattern("{sigfile}")),
    output: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    wildcard_constraints:
        input_type="genomic|genomic-reps"
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))


localrules: protein_signames_to_file 
rule protein_signames_to_file:
    input: 
        csv=os.path.join(out_dir, ".{basename}.{input_type}.check_csv"),
        chpt=os.path.join(out_dir, ".{basename}.{input_type}.check_proteins"),
        sigs=ancient(Check_Protein_Sigs("{sigfile}")),
    output: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    wildcard_constraints:
        input_type="protein|protein-reps"
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                outF.write(full_filename + "\n")

localrules: sigs_to_zipfile
rule sigs_to_zipfile:
    input: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    output: os.path.join(out_dir, "databases", "{basename}.{input_type}.zip")
    params:
        compression=config.get("gzip_compression", 9)
    log: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.sigs-to-zipfile.log")
    benchmark: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.sigs-to-zipfile.benchmark")
    conda: "envs/sourmash4.yml"
    shell:
        """
        python sigs-to-zipfile.py --compression {params.compression} --sig-pathlist {input} {output} 2> {log}
        """

localrules: sigs_to_ksize_zipfile
rule sigs_to_ksize_zipfile:
    input: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    output: os.path.join(out_dir, "ksize_databases", "{basename}.{input_type}.k{ksize}.zip")
    params:
        compression=config.get("gzip_compression", 9)
    log: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.k{ksize}.sigs-to-zipfile.log")
    benchmark: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.k{ksize}.sigs-to-zipfile.benchmark")
    conda: "envs/sourmash4.yml"
    shell:
        """
        python sigs-to-zipfile.py --compression {params.compression} --ksize {wildcards.ksize} --sig-pathlist {input} {output} 2> {log}
        """

