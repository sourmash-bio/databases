import os, sys
import pandas as pd
import glob

configfile: "defaults.yml"
#configfile: "index_config.yml" # original gtdb config

out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")
data_dir = config['data_dir'].rstrip('/')
basename = config["basename"]

def sanitize_path(path):
    # expand `~`, get absolute path
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    return path

## this file must have at least two columns: accession,filename
def read_samples(samples_file, data_dir):
    samples = pd.read_csv(samples_file, dtype=str, sep=",", header=0)
    samples["accession"] = samples["accession"].str.split(" ", expand=True)[0]
    if "signame" not in samples.columns:
        if "species" in samples.columns:
            samples['signame'] = samples["accession"] + " " + samples["species"]
        else:
            samples['signame'] = samples["accession"]
    samples.set_index("accession", inplace=True)
    
    # Now, verify that all genome files exist
    data_dir = sanitize_path(data_dir)
    sample_list = samples["filename"].tolist()
    for filename in sample_list:
        fullpath = os.path.join(data_dir, filename)
        if not os.path.exists(fullpath):
            print(f'** ERROR: genome file {filename} does not exist in {data_dir}')
    return samples

# figure out inputs
dna_input = config.get("dna_input", True)
protein_input = config.get("protein_input", False)
singleton = config.get("singleton", False)
alphabet_info = config["alphabet_info"]
sample_names = []
input_types=[]

if dna_input:
    genome_info = read_samples(config["genomes_csv"], data_dir)
    sample_names = genome_info.index.tolist()
    input_types +=["dna-input"]
if protein_input:
    protein_info = read_samples(config["proteins_csv"], data_dir)
    input_types +=["protein-input"]
    if not sample_names:
        sample_names = protein_info.index.tolist()
elif any(alpha in ["protein", "dayhoff", "hp"] for alpha in alphabet_info):
    print("Error: protein alphabets found in the desired sbt alphabets. Please set 'protein_input: True' in your config and provide a csv with protein input ('proteins_csv').")
    sys.exit(-1)

onstart:
    print("------------------------------")
    print("   Index from lineages csv")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

index_info= []

wildcard_constraints:
    alphabet="\w+",
    ksize="\d+"

for alphabet, info in alphabet_info.items():
    aks = expand("{alpha}-k{ksize}-scaled{scaled}", alpha=alphabet, ksize=info["ksizes"], scaled=info["scaled"])
    index_info.extend(aks)

index_types = ["sbt.zip"]
make_lca = config.get("index_lca", False)
if make_lca:
    index_types.append("lca.json.gz")


rule all:
    input: 
        expand(os.path.join(out_dir, "{input_type}", "{basename}.signatures.txt"), basename=basename, input_type = input_types),
        expand(os.path.join(out_dir, "index", "{basename}.{idx_info}.{idx_type}"), basename=basename, idx_info=index_info, idx_type=index_types),


## sketching rules ##
def build_sketch_params(output_type):
    sketch_cmd = ""
    if output_type == "nucleotide":
        ksizes = config["alphabet_info"]["nucleotide"].get("ksizes", config["alphabet_defaults"]["nucleotide"]["ksizes"])
        # just build signatures at the highest resolution (lowest scaled value). Index can downsample for us.
        scaled = min(config["alphabet_info"]["nucleotide"].get("scaled", config["alphabet_defaults"]["nucleotide"]["scaled"]))
        # maybe don't track abundance?
        sketch_cmd = "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
        return sketch_cmd
    elif output_type == "protein":
        for alpha in ["protein", "dayhoff", "hp"]:
            if alpha in config["alphabet_info"].keys():
                ## if ksizes aren't given, sketch protein, dayhoff, hp at the ksizes from default config
                ksizes = config["alphabet_info"][alpha].get("ksizes", config["alphabet_defaults"][alpha]["ksizes"])
                scaled = config["alphabet_info"][alpha].get("scaled", config["alphabet_defaults"][alpha]["scaled"])
            else:
                ksizes = config["alphabet_defaults"][alpha]["ksizes"]
                scaled = config["alphabet_defaults"][alpha]["scaled"]
            sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

def build_signame_cmd(sample, input_type):
    if not singleton:
        if input_type == "nucleotide":
            signame = genome_info.at[sample, 'signame']
        elif input_type == "protein":
            sgname = protein_info.at[sample, 'signame'],
        signame_cmd =  f" --name {signame:q}"
        return signame_cmd
    else:
        return ""

rule sourmash_sketch_nucleotide_input:
    input: 
        lambda w: os.path.join(data_dir, genome_info.at[w.sample, 'filename'])
    output:
        os.path.join(out_dir, "dna-input", "signatures", "{sample}.sig"),
    params:
        sketch_params = build_sketch_params("nucleotide"),
        signame_cmd = lambda w: build_signame_cmd(w.sample, "nucleotide"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    group: "sigs"
    log: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.log")
    benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} {params.signame_cmd} {input}  2> {log}
        """
    
if protein_input:
    rule sourmash_sketch_protein_input:
        input: lambda w: os.path.join(data_dir, protein_info.at[w.sample, 'filename'])
        output:
            os.path.join(out_dir, "protein-input", "signatures", "{sample}.sig"),
        params:
            sketch_params = build_sketch_params("protein"),
            signame_cmd = lambda w: build_signame_cmd(w.sample, "protein"),
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *1000,
            runtime=1200,
        group: "sigs"
        log: os.path. join(logs_dir, "sourmash_sketch_prot_input", "{sample}.sketch.log")
        benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_prot_input", "{sample}.sketch.benchmark")
        conda: "envs/sourmash-dev.yml"
        shell:
            """
            sourmash sketch protein {params.sketch_params} -o {output} {params.signame_cmd} {input} 2> {log}
            """

localrules: signames_to_file

checkpoint signames_to_file:
    input:  expand(os.path.join(out_dir, "{{input_type}}", "signatures", "{sample}.sig"), sample=sample_names),
    output: os.path.join(out_dir, "{input_type}", "{basename}.signatures.txt")
    group: "sigs"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


def get_siglist(w):
    if w.alphabet in ["protein", "dayhoff", "hp"]:
        return os.path.join(out_dir, "protein-input", "{basename}.signatures.txt")
    else:
        return os.path.join(out_dir, "dna-input", "{basename}.signatures.txt")


rule index_sbt:
    input: get_siglist 
    output: os.path.join(out_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.sbt.zip"),
    threads: 1
    params:
        alpha_cmd = lambda w: config["alphabet_defaults"][w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(config["alphabet_defaults"][w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-sbt.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-sbt.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash index {output} --ksize {params.ksize} \
        --scaled {wildcards.scaled} {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

def get_lineages(w):
    if w.alphabet in ["protein", "dayhoff", "hp"]:
        return config["proteins_csv"]
    else:
        return config["genomes_csv"]

rule index_lca:
    input:
        sigfile=get_siglist,
        lineages=get_lineages
    output: os.path.join(out_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.lca.json.gz"),
    threads: 1
    params:
        alpha_cmd = lambda w: config["alphabet_defaults"][w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(config["alphabet_defaults"][w.alphabet]["ksize_multiplier"]),
        report = lambda w: os.path.join(out_dir,"index", f"{w.basename}.{w.alphabet}-k{w.ksize}-scaled{w.scaled}.lca.report"),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    log: os.path.join(logs_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-lca.log")
    benchmark: os.path.join(benchmarks_dir, "index", "{basename}.{alphabet}-k{ksize}-scaled{scaled}.index-lca.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash lca index --ksize {params.ksize} \
          --scaled {wildcards.scaled} --from-file {input.sigfile} \
          --start-column 3 --require-taxonomy --report {params.report} \
          {params.alpha_cmd} {input.lineages} {output} 2> {log} \
        """
