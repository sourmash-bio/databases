DOMAINS=['fungi', 'viral']
GENBANK_INPUTS = {}
for domain in DOMAINS:
    GENBANK_INPUTS[domain] = [l for l in shell('find /home/irber/ncbi/genbank/{domain} -iname "*_genomic.fna.gz"', iterable=True) if l]

import os.path
def relpath(x):
    return os.path.relpath(x, '/home/irber/ncbi')

#with open('genbank_inputs', 'r') as f:
#    GENBANK_INPUTS = [l.strip() for l in f.readlines()]

rule all:
    input:
        expand("outputs/trees/scaled/{db}-{domain}-d{nchildren}-k{ksize}.sbt.json", db=['genbank'], nchildren=[10], ksize=[21], domain=DOMAINS)

rule tetramer_sigs:
    input: "{db}/{group}/{id}/{filename}_genomic.fna.gz"
    output: "outputs/sigs/4.5.mers/{db}/{group}/{id}/{filename}_genomic.fna.gz.sig"
    params: filename="{filename}"
    threads: 1
	shell: """
		mkdir -p `dirname {output}`
		sourmash compute -k 4,5 \
                         -n 2000 \
                         --track-abundance \
                         --name-from-first \
                         -o {output} \
                         {input}
	"""

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

def sbt_inputs(w):
    if w.db == 'refseq':
        return expand("outputs/sigs/{config}/{ids}.sig", config=w.config, ids=REFSEQ_INPUTS)
    elif w.db == 'genbank':
        x = expand("outputs/sigs/{config}/{ids}.sig", config=w.config, ids=[ relpath(x) for x in GENBANK_INPUTS[w.domain] ])
	print(x)
	return x

    print(repr(w), w.db)


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
