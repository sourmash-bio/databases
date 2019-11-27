import os.path
import random

configfile: "config.yml"

# number of signatures to calculate in a single run
BATCHSIZE=5000

# load in list of genomes by domain
DOMAINS=config['domains']

# load in files for each domain
GENBANK_INPUTS = {}
for domain in DOMAINS:
    GENBANK_INPUTS[domain] = [ x.strip() for x in open('domain-{}.txt'.format(domain)) ]

# function to remove prefix
def relpath(x):
    return os.path.relpath(x, config['genome_location'])

def list_all_sigs():
    "Get a list of signatures that may need to be calculated."
    all_sigs = []

    if os.path.exists('sigs-todo.txt'):
        all_sigs = [ x.strip() for x in open('sigs-todo.txt') ]
        print('loaded {} sigs from sigs-todo.txt'.format(len(all_sigs)))
    else:
        print('loading sigs from individual domains...')
        for domain in DOMAINS:
            ids = [ relpath(x) for x in GENBANK_INPUTS[domain] ]
            all_sigs.extend(expand("outputs/sigs/scaled/{ids}.sig", ids=ids))
    print('calculating up to {} sigs'.format(len(all_sigs)))

    random.shuffle(all_sigs)

    return all_sigs

def calc_undone_sigs():
    "Collect up to BATCHSIZE undone signatures."
    global full_list_sigs
    global done_sigs

    batch = []
    n = 0
    for n, sigfile in enumerate(full_list_sigs):
        if os.path.exists(sigfile):
            done_sigs.add(sigfile)
        else:
            batch.append(sigfile)
            
        if len(batch) >= BATCHSIZE:
            break

    print('looked at {}, calculating {}'.format(n, len(batch)))

    if done_sigs:
        remaining = set(full_list_sigs) - done_sigs
        print('removed known done sigs, from {} to {}'.format(len(full_list_sigs), len(remaining)))
        with open('sigs-todo.txt', 'wt') as fp:
            fp.write("\n".join(remaining))

    print(batch[:10])

    return batch

###

# all output sigs
full_list_sigs = list_all_sigs()
done_sigs = set()

# rule 'all' builds specific SBTs.
rule all:
    input:
        expand("outputs/trees/scaled/{db}-{domain}-d{nchildren}-x{bfsize}-k{ksize}.sbt.json", db=['genbank'], nchildren=[2, 10], ksize=config['db_ksizes'], domain=DOMAINS, bfsize=["1e4", "1e5", "1e6"]),
        expand("outputs/lca/scaled/{db}-{domain}-k{ksize}-scaled10k.lca.json.gz", db=['genbank'], domain=DOMAINS, ksize=config['db_ksizes'])

# build scaled signatures needed for the SBTs.
rule scaled_sigs:
    input: "/home/irber/ncbi/{db}/{group}/{id}/{filename}_genomic.fna.gz"
    output: "outputs/sigs/scaled/{db}/{group}/{id}/{filename}_genomic.fna.gz.sig"
    params: filename="{filename}"
    threads: 1
    shell: """
		mkdir -p `dirname {output}`
		sourmash compute -k 21,31,51 \
                         --scaled 2000 \
                         --track-abundance \
                         --name-from-first \
                         -o {output} \
                         {input}
	"""

# return all of the inputs needed for a specific SBT.
def sbt_inputs(w):
    if w.db == 'refseq':
        return expand("outputs/sigs/{config}/{ids}.sig", config=w.config, ids=REFSEQ_INPUTS)
    elif w.db == 'genbank':
        x = expand("outputs/sigs/{config}/{ids}.sig", config=w.config, ids=[ relpath(x) for x in GENBANK_INPUTS[w.domain] ])
	print('gathered {} sbt inputs for {} / {}'.format(len(x), w.db, w.domain))
	return x

    print(repr(w), w.db)

# rule 'all_sigs' will build all specific signature files.
# (maybe this will run faster than relying on SBT rule to specify?)

rule all_sigs:
    input:
        calc_undone_sigs()

# build actual SBT!
rule sbt_tree:
    input: sbt_inputs
    output: "outputs/trees/{config}/{db}-{domain}-d{nchildren}-x{bfsize}-k{ksize}.sbt.json"
    params:
        ksize="{ksize}",
        db="{db}",
        config="{config}",
        nchildren="{nchildren}",
	domain="{domain}",
        bfsize="{bfsize}"
    shell: """
		mkdir -p `dirname {output}`
        sourmash index -k {params.ksize} \
                           -d {params.nchildren} \
                           -x {params.bfsize} \
                           --traverse-directory \
                           {output} outputs/sigs/{params.config}/{params.db}/{params.domain}
    """

# build LCA
rule lca_db:
    input:
        sbt_inputs,
        "domain-{domain}.acc.lineages.csv"
    output:
        "outputs/lca/{config}/{db}-{domain}-k{ksize}-scaled10k.lca.json.gz",
    params:
        ksize="{ksize}",
        db="{db}",
	domain="{domain}",
        config="{config}"
    shell: """
		mkdir -p `dirname {output}`
        sourmash lca index -k {params.ksize} \
                           --scaled 10000 \
                           --traverse-directory -C 3 --split-identifiers \
                           domain-{params.domain}.acc.lineages.csv \
                           {output} outputs/sigs/{params.config}/{params.db}/{params.domain}
    """

# deprecated
rule upload_all_to_hpcc:
    input:
        expand("outputs/tar_trees/scaled/{db}-k{ksize}.tar.gz", db=['refseq', 'genbank'], nchildren=[2], ksize=[21, 31, 51]),
        expand("outputs/tar_trees/4.5.mers/{db}-k{ksize}.tar.gz", db=['refseq', 'genbank'], nchildren=[2], ksize=[4, 5])

rule upload_one_to_hpcc:
    input: "outputs/trees/{config}/{db}-k{ksize}.sbt.json"
    output: "outputs/tar_trees/{config}/{db}-k{ksize}.tar.gz"
    params:
        ksize="{ksize}",
        db="{db}",
        config="{config}",
    shell: """
        mkdir -p `dirname {output}`
        export CUR_DIR=`pwd`
        cd `dirname {input}`
        tar zcf ${{CUR_DIR}}/{output} \
            {params.db}-k{params.ksize}.sbt.json \
            .sbt.{params.db}-k{params.ksize}/
        cd -
        scp {output} hpcc:/mnt/research/ged/irberlui/sbts/
    """
