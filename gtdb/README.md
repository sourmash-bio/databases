# build gtdb databases for sourmash

## First, test your installation etc.

0. Install `sourmash` v2 or v3.

1. Download the [GTDB-Tk reference data](https://github.com/Ecogenomics/GtdbTk#gtdb-tk-reference-data) (currently release89)

2. Edit `config.yml` in this directory to point to `release89/taxonomy/gtdb_taxonomy.tsv`.

3. Make a directory that contains a ~dozen of the genomes from `release89/fastani/datbase`, and set `genomes_location` in config.yml to that directory.

4. Set `sigs_output_location` to a local directory (you can leave it as-is if you like).

5. Run `snakemake`. You can use `-j 8` if you like.

At the end you should have one `*.sig` file in your `sigs_output_location` for each genome in `genomes_location`, and six output files in the current directory:

```
gtdb-release89-k21.lca.json.gz
gtdb-release89-k31.lca.json.gz
gtdb-release89-k51.lca.json.gz
gtdb-release89-k21.sbt.json
gtdb-release89-k31.sbt.json
gtdb-release89-k51.sbt.json
```

The first three are standalone LCA databases for sourmash, and the
last three are SBT databases for sourmash; see also the
`.sbt.*` hidden directories.

## Build the full databases

1. Copy `config.yml` to `config-release89.yml`

2. Comment out `prefix` with a `#` at the beginning of the line.

3. Change `genomes_location` to point to `release89/fastani/database`.

4. Change `sigs_output_location` to point to another directory.

5. Run `snakemake`.

With `snakemake -j 8`, building the 25,000 signatures takes about 3 hours
on the farm cluster at UC Davis, and each database takes about an hour to
build after that. Memory requirements are far under 100 GB total.
