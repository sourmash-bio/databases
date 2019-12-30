~/2019-sourmash-gtdb/find-oddities.py gtdb-release89-k51.lca.json.gz --lowest-rank=superkingdom --minimum-hashes=10 --prefix=oddities-k51 > oddities-k51.txt
~/2019-sourmash-gtdb/find-oddities-examine.py oddities-k51.csv ~/gtdbtk/release89/fastani/database --percent-threshold=95 --length-threshold=0 > oddities-k51.examine.txt

~/2019-sourmash-gtdb/find-oddities.py gtdb-release89-k21.lca.json.gz --lowest-rank=superkingdom --minimum-hashes=10 --prefix=oddities-k21 > oddities-k21.txt
~/2019-sourmash-gtdb/find-oddities-examine.py oddities-k21.csv ~/gtdbtk/release89/fastani/database --percent-threshold=95 --length-threshold=0  > oddities-k21.examine.txt
