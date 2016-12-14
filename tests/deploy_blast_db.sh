mkdir -p local/blast
cat raw/*.fna > local/blast/bacteria.fa
makeblastdb -dbtype nucl -in local/blast/bacteria.fa
