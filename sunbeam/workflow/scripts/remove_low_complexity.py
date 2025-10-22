import traceback
from typing import TextIO


def f(log: TextIO):
    import os
    import subprocess as sp

    input_reads = snakemake.input.reads  # type: ignore
    output_reads = snakemake.output.reads  # type: ignore
    output_report = snakemake.output.report  # type: ignore
    kmer_size = snakemake.params.kmer_size  # type: ignore
    min_kscore = snakemake.params.min_kscore  # type: ignore
    threads = snakemake.threads  # type: ignore
    compression = snakemake.params.compression  # type: ignore

    args = [
        "heyfastq",
        "filter-kscore",
        "--input",
        *input_reads,
        "--output",
        *output_reads,
        "--report",
        output_report,
        "--kmer-size",
        str(kmer_size),
        "--min-kscore",
        str(min_kscore),
        "--threads",
        str(threads),
    ]

    os.environ["HFQ_GZIP_COMPRESSION"] = str(compression)
    res = sp.run(args, capture_output=True, text=True)
    log.write(res.stdout)
    log.write(res.stderr)
    if res.returncode != 0:
        raise RuntimeError(
            f"heyfastq filter-kscore failed with exit code {res.returncode}"
        )


log_f = snakemake.log[0]  # type: ignore
with open(log_f, "w") as log:
    log.write("Initialized log and starting script...\n")
    stderr_capture = None
    try:
        f(log)
    except BaseException as e:
        log.write(f"Error during run: {e}\n")
        log.write(traceback.format_exc())
        raise
    else:
        log.write("Completed successfully.\n")
