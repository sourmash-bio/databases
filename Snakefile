# load in list of genomes by domain
DOMAINS=['fungi', 'viral', 'bacteria', 'archaea']
GENBANK_INPUTS = {}
for domain in DOMAINS:
#    GENBANK_INPUTS[domain] = [l for l in shell('find /home/irber/ncbi/genbank/{domain} -iname "*_genomic.fna.gz"', iterable=True) if l]
    GENBANK_INPUTS[domain] = [ x.strip() for x in open('domain-{}.txt'.format(domain)) ]

# function to remove prefix
import os.path
def relpath(x):
    return os.path.relpath(x, '/home/irber/ncbi')

def list_all_sigs():
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

    import random
    random.shuffle(all_sigs)

    return all_sigs

def calc_undone_sigs():
    global full_list_sigs
    global done_sigs

    BATCHSIZE=5000
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

# calculate all output sigs / hardcoded for now
full_list_sigs = list_all_sigs()
done_sigs = set()

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

rule all_sigs:
    input:
        calc_undone_sigs()

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
