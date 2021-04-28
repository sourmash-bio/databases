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

siglist_types = []
if config.get("build_genomic", False):
    siglist_types.append("genomic")
    if config.get("build_representative_sets", False):
        siglist_types.append("genomic-reps")
if config.get("build_protein", False):
    siglist_types.append("protein")
    if config.get("build_representative_sets", False):
        siglist_types.append("protein-reps")

# set signature directories
existing_genomic_sigdir = config["genomic"]["sig_dir"]
new_genomic_sigdir = os.path.join(out_dir, "genomic", "signatures")

existing_protein_sigdir = config["protein"]["sig_dir"]
new_protein_sigdir = os.path.join(out_dir, "protein", "signatures")

# check params are in the right format
for alpha, info in config["alphabet_info"].items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        config["alphabet_info"][alpha]["scaled"] = [scaled]
    if not isinstance(ksize, list):
        config["alphabet_info"][alpha]["ksize"] = [ksize]

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

    def find_sigs(self, genomic=False, protein=False):
        global tax_info
        global sketch_genomic
        global sketch_protein
        sketch_genomic = pd.DataFrame()
        sketch_protein = pd.DataFrame()

        tax_info = pd.read_csv(f"{out_dir}/{basename}.taxonomy.csv", header=0)
        #For metadata info: DtypeWarning: Columns (61,65,74,82,83) have mixed types.Specify dtype option on import or set low_memory=False
        metadata_info = pd.read_csv(f"{out_dir}/{basename}.metadata.csv.gz", header=0, low_memory=False) 
        representative_accessions = metadata_info[metadata_info["gtdb_representative"] == "t"]["accession"].str.replace("GB_", "").str.replace("RS_", "")
        
        # build dictionaries of accessions we need to download and sketch (genomic, protein)
        if genomic:
            print("finding genomic sigs...")
            # build sigfile basename
            tax_info["genomic_sigfile"] = tax_info["ident"] + ".sig"
            # find all sigpaths:
            tax_info["genomic_sigfile"] = tax_info["genomic_sigfile"].apply(update_sigpath, sigdir=existing_genomic_sigdir, new_sigdir=new_genomic_sigdir)
            
            # DF for just the files that need to be downloaded and sketched
            sketch_genomic = tax_info.loc[~tax_info["genomic_sigfile"].apply(sigfile_exists)]
            sketch_genomic["signame"] = sketch_genomic["ident"] + " " + sketch_genomic["species"].str.replace("s__", "")
            
            # find representative genome files
            representative_genomic = tax_info.loc[tax_info["ident"].isin(representative_accessions)]
            
            sketch_genomic.set_index("ident", inplace=True)
            print(f"found {len(sketch_genomic)} genomic accessions that need to be downloaded and sketched") 
            
            all_genomic_sigfiles = tax_info["genomic_sigfile"]
            rep_genomic_sigfiles = representative_genomic["genomic_sigfile"]
            
            return all_genomic_sigfiles, rep_genomic_sigfiles
            

        if protein:
            print("finding protein sigs...")
            tax_info["protein_sigfile"] = tax_info["ident"] + ".sig"
            # find all sigpaths:
            tax_info["protein_sigfile"] = tax_info["protein_sigfile"].apply(update_sigpath, sigdir=existing_protein_sigdir, new_sigdir=new_protein_sigdir)
            
            # DF for just the files that need to be downloaded and sketched
            sketch_protein = tax_info[~tax_info["protein_sigfile"].apply(sigfile_exists)]
            sketch_protein["signame"] = sketch_protein["ident"] + " " + sketch_protein["species"].str.replace("s__", "")
            print(f"found {len(sketch_protein)} protein accessions that need to be downloaded and sketched") 
            sketch_protein.set_index("ident", inplace=True)

            # ONLY use genome accessions with protein assemblies
            all_protein_sigfiles = tax_info[tax_info["ident"].isin(sketch_protein.index)]["protein_sigfile"]
            #all_protein_sigfiles = tax_info["protein_sigfile"]
            
            # find representative protein files
            representative_protein = tax_info.loc[tax_info["ident"].isin(representative_accessions)]
            rep_protein_sigfiles = representative_protein["protein_sigfile"]
            
            return all_protein_sigfiles, rep_protein_sigfiles
        

    def __call__(self, w):
        global checkpoints
        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.

        checkpoints.check_csv.get(**w)
        
        # this works, but there is probably a better way to do it...
        genomic_sigfiles, protein_sigfiles, rep_genomic_sigfiles, rep_protein_sigfiles=[],[],[],[]
        if config["build_genomic"]:
            genomic_sigfiles, rep_genomic_sigfiles = self.find_sigs(genomic=True)
        if config["build_protein"]:
            protein_sigfiles, rep_protein_sigfiles = self.find_sigs(protein=True)
        
        sigfiles = {"genomic": genomic_sigfiles, 
                    "protein": protein_sigfiles,
                    "genomic-reps": rep_genomic_sigfiles, 
                    "protein-reps": rep_protein_sigfiles }

        pattern = expand(self.pattern, sigfile=sigfiles[w.input_type], **w)
        return pattern


rule all:
    input: expand(os.path.join(out_dir, "{basename}.{input_type}.zip"), basename=basename, input_type=siglist_types)

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
        taxonomy=os.path.join(out_dir, f"{basename}.taxonomy.csv"),
        metadata=os.path.join(out_dir, f"{basename}.metadata.csv.gz")
    output: touch(f"{out_dir}/.check_csv")

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
        # TODO: might check here if it was a 404 or 5xx, assuming 404
        raise ValueError(f"Can't find URL for {accession}, tried {url}")

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
        signame = lambda w: sketch_genomic.at[w.accession, "signame"],
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
        signame = lambda w: sketch_protein.at[w.accession, "signame"],
        url_path = lambda w: url_for_accession(w.accession, protein=True)
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
        """

localrules: signames_to_file
rule signames_to_file:
    input: 
        csv=f"{out_dir}/.check_csv",
        sigs=ancient(Checkpoint_MakePattern("{sigfile}")),
    output: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                outF.write(full_filename + "\n")

localrules: sigs_to_zipfile
rule sigs_to_zipfile:
    input: os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    output: os.path.join(out_dir, "{basename}.{input_type}.zip")
    params:
        compression=config.get("gzip_compression", 9)
    log: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.sigs-to-zipfile.log")
    benchmark: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.{input_type}.sigs-to-zipfile.benchmark")
    conda: "envs/sourmash4.yml"
    shell:
        """
        python sigs-to-zipfile.py --compression {params.compression} --sig-pathlist {input} {output} 2> {log}
        """

