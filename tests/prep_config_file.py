import sys
import yaml
import argparse

def main():

    parser = argparse.ArgumentParser(description="Modifies test config with testing values")
    parser.add_argument(
        "--config", help="Test config file", type=argparse.FileType("r"),
        default=sys.stdin)

    args = parser.parse_args()
    config = yaml.load(args.config)

    config['all']['filename_fmt'] = "PCMP_{sample}_{rp}.fastq"
    config['qc']['human_index_fp'] = "indexes/human.fasta"
    config['qc']['phix_index_fp'] = "indexes/phix174.fasta"
    config['classify']['kraken_db_fp'] = "mindb"
    config['assembly']['cap3_fp'] = "local/CAP3"
    config['blastdbs']['root_fp'] = "local/blast"
    config['blastdbs']['nucleotide']['bacteria'] = 'bacteria.fa'
    config['mapping']['genomes_fp'] = "indexes"
    config['mapping']['igv_fp'] = "local/IGV/igv.sh"

    sys.stdout.write(yaml.dump(config))

if __name__ == "__main__":
    main()
