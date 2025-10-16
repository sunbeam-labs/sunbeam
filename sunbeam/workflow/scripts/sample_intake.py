import os


log_fp = snakemake.log[0]  # type: ignore
with open(log_fp, "w") as log:
    try:
        log.write("Starting ...\n")

        input_reads = snakemake.input[0]  # type: ignore
        output_reads = snakemake.output[0]  # type: ignore

        assert input_reads.endswith(".fastq.gz")
        os.symlink(input_reads, output_reads)
    except BaseException as e:
        log.write(f"Error during run: {e}\n")
        raise e
