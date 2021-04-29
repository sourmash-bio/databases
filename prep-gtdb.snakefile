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

    def find_sigs(self, genomic=False, protein=False):
        global tax_info
        global sketch_genomic
        global sketch_protein
        global representative_accessions
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

            all_protein_sigfiles = tax_info["protein_sigfile"]
            #all_protein_sigfiles = tax_info[tax_info["ident"].isin(sketch_protein.index)]["protein_sigfile"]
            
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


class Check_Protein_Sigs:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_accessions_and_update_sigpaths(self):
        with open(f"{data_dir}/{basename}.protein-reps.prodigal-siglist.txt", "rt") as fp:
            accessions = [ x.rstrip().rsplit(".sig") for x in fp ]
            for accession in accessions:
                taxinfo.at[accession, "protein_sigfile"] = f"{prodigal_sigdir}/{accession}.sig"
            
            all_protein_sigfiles = tax_info["protein_sigfile"]
            
            # find representative protein files
            representative_protein = tax_info.loc[tax_info["ident"].isin(representative_accessions)]
            rep_protein_sigfiles = representative_protein["protein_sigfile"]
            
        return all_protein_sigfiles, rep_protein_sigfiles

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_proteins'; this will trigger an
        # exception until that rule has been run.
        #checkpoints.check_proteins.get(**w)
        checkpoints.check_protein_sigfiles_for_empties.get(**w)
        protein_sigfiles, rep_protein_sigfiles = self.get_accessions_and_update_sigpaths()

        sigfiles = {"protein": protein_sigfiles,
                    "protein-reps": rep_protein_sigfiles}

        pattern = expand(self.pattern, sigfile=sigfiles[w.input_type], **w)
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
        #raise ValueError(f"Can't find URL for {accession}, tried {url}")
        return ""

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
        url_path = lambda w: url_for_accession(w.accession, protein=True),
        #error_file = os.path.join(new_protein_sigdir, "{accession}.sig.err"),
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
    # could touch output file in case it fails
    #if [[ -s {output} ]] || touch {params.error_file}

checkpoint check_protein_sigfiles_for_empties:
    input:
        csv=f"{out_dir}/.check_csv",
        sigs=ancient(Checkpoint_MakePattern("{sigfile}")),
    output: os.path.join(out_dir, "{basename}.{input_type}.prodigal-siglist.txt")
    wildcard_constraints:
        input_type="protein|protein-reps"
    run:
        with open(output, 'w') as out:
            for sig in input.sigs:
                if os.path.getsize(str(sig)) == 0:
                    out.write(f"{str(sig)}\n")
                

# make the checkpoint work
#localrules: check_proteins
#checkpoint check_proteins:
#    input:
#        expand(os.path.join(out_dir, "{basename}.{input_type}.prodigal-siglist.txt"), basename = basename, input_type="protein-reps")
#        #os.path.join(data_dir, f"{basename}.prodigal-protein.txt")
#    wildcard_constraints:
#        input_type="protein|protein-reps"
#    output: touch(f"{out_dir}/.check_proteins")

rule download_genomes_for_failed_protein_sigs:
    output: temp(os.path.join(data_dir, "{accession}_genomic.fna")) # output file marked as temp is deleted after all rules that use it as an input are completed
    params:
        url_path = lambda w: url_for_accession(w.accession),
    shell:
        """
        curl -s {params.url_path} | zcat > {output} 2> {log}
        """

rule prodigal_genomes_for_failed_protein_sigs:
    input: os.path.join(data_dir, "{accession}_genomic.fna")
    output: temp(os.path.join(data_dir, "{accession}.prodigal.faa"))
    conda: "envs/prodigal-env.yml"
    log: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.log")
    benchmark: os.path.join(logs_dir, "prodigal", "{accession}.prodigal.benchmark")
    shell:
        """
         prodigal -i {input} -a {output.proteins} 2> {log}
        """

rule sketch_prodigal_genomes_for_failed_protein_downloads:
    input: os.path.join(data_dir, "{accession}.prodigal.faa")
    output: os.path.join(prodigal_sigdir,  "{accession}.prodigal.sig"),
    params:
        sketch_params = make_protein_param_str(),
        signame = lambda w: sketch_protein.at[w.accession, "signame"],
        url_path = lambda w: url_for_accession(w.accession, protein=True),
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
        csv=f"{out_dir}/.check_csv",
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
        #check_proteins=f"{out_dir}/.check_proteins",
        check_proteins=os.path.join(out_dir, "{basename}.{input_type}.prodigal-siglist.txt"),
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

