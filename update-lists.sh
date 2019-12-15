for domain in fungi viral bacteria archaea; do find /home/irber/ncbi/genbank/$domain -iname "*_genomic.fna.gz" > domain-$domain.txt; done
