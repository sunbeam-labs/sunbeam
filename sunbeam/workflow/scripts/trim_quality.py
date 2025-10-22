import traceback
from typing import TextIO


def f(log: TextIO):
    import os
    import subprocess as sp

    input_reads = snakemake.input.reads  # type: ignore
    output_reads = snakemake.output.reads  # type: ignore
    output_report = snakemake.output.report  # type: ignore
    window_width, window_threshold = snakemake.params.window  # type: ignore
    start_threshold = snakemake.params.start_threshold  # type: ignore
    end_threshold = snakemake.params.end_threshold  # type: ignore
    min_length = snakemake.params.min_length  # type: ignore
    threads = snakemake.threads  # type: ignore
    compression = snakemake.params.compression  # type: ignore

    args = [
        "heyfastq",
        "trim-qual",
        "--input",
        *input_reads,
        "--output",
        *output_reads,
        "--report",
        output_report,
        "--window-width",
        str(window_width),
        "--window-threshold",
        str(window_threshold),
        "--start-threshold",
        str(start_threshold),
        "--end-threshold",
        str(end_threshold),
        "--min-length",
        str(min_length),
        "--threads",
        str(threads),
    ]

    os.environ["HFQ_GZIP_COMPRESSION"] = str(compression)
    res = sp.run(args, capture_output=True, text=True)
    log.write(res.stdout)
    log.write(res.stderr)
    if res.returncode != 0:
        raise RuntimeError(f"heyfastq trim-qual failed with exit code {res.returncode}")


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
