# load in list of genomes by domain
DOMAINS=['fungi', 'viral'] #, 'bacteria', 'archaea']
GENBANK_INPUTS = {}
for domain in DOMAINS:
#    GENBANK_INPUTS[domain] = [l for l in shell('find /home/irber/ncbi/genbank/{domain} -iname "*_genomic.fna.gz"', iterable=True) if l]
    GENBANK_INPUTS[domain] = [ x.strip() for x in open('domain-{}.txt'.format(domain)) ]

# function to remove prefix
import os.path
def relpath(x):
    return os.path.relpath(x, '/home/irber/ncbi')

# rule 'all' builds specific SBTs.
rule all:
    input:
        expand("outputs/trees/scaled/{db}-{domain}-d{nchildren}-k{ksize}.sbt.json", db=['genbank'], nchildren=[10], ksize=[21], domain=DOMAINS)

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
	print(x)
	return x

    print(repr(w), w.db)

# rule 'all_sigs' will build all specific signature files.
# (maybe this will run faster than relying on SBT rule to specify?)

def list_all_sigs():
    all_sigs = []
    for domain in DOMAINS:
        ids = [ relpath(x) for x in GENBANK_INPUTS[domain] ]
        all_sigs.extend(expand("outputs/sigs/scaled/{ids}.sig", ids=ids))
    return all_sigs

rule all_sigs:
    input:
        list_all_sigs()

# build actual SBT!
rule sbt_tree:
    input: sbt_inputs
    output: "outputs/trees/{config}/{db}-{domain}-d{nchildren}-k{ksize}.sbt.json"
    threads: 32
    params:
        ksize="{ksize}",
        db="{db}",
        config="{config}",
        nchildren="{nchildren}",
    shell: """
		mkdir -p `dirname {output}`
        sourmash index -k {params.ksize} \
                           -d {params.nchildren} \
                           -x 1e6 \
                           --traverse-directory \
                           {output} outputs/sigs/{params.config}/{params.db}
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
