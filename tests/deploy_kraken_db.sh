set -e
set -x

KRAKEN_DB_NAME="/home/chunyu/sunbeam/tests/local/db/bacteria"

# download the taxonomy data
kraken-build --db $KRAKEN_DB_NAME --download-taxonomy

# download the genome_urls.txt
snakemake -s snakefile_kraken_db

# download the genome sequences
snakemake -s snakefile_kraken_db download_group

# add to kraken db
snakemake -s snakefile_kraken_db add_group_to_kraken_db

# build the database
kraken-build --build --db $KRAKEN_DB_NAME --threads 8 --jellyfish-hash-size 3200M

# clean the database
kraken-build --clean --db $KRAKEN_DB_NAME

echo $"Finish building kraken database"



