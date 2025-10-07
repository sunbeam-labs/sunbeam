import subprocess as sp
import traceback
from pathlib import Path
from typing import Optional, TextIO


def f(log: TextIO) -> Optional[sp.CompletedProcess[str]]:
    input_reads = snakemake.input.reads  # type: ignore
    output_reads = snakemake.output.reads  # type: ignore
    failed_reads = snakemake.output.fail  # type: ignore
    json = snakemake.output.json  # type: ignore
    html = snakemake.output.html  # type: ignore
    adapter = snakemake.params.adapter  # type: ignore
    threads = snakemake.threads  # type: ignore

    adapter_fp = Path(adapter)
    assert (
        adapter_fp.is_file()
    ), f"Adapter file {adapter} does not exist. Point config to an empty file to skip this step."
    if adapter_fp.stat().st_size == 0:
        log.write(
            f"Warning: Adapter file {adapter} is empty. No adapter removal will be performed.\n"
        )
        return None

    args = [
        "fastp",
        "--thread",
        str(threads),
        "-i",
        str(input_reads[0]),
        "-o",
        str(output_reads[0]),
        "--detect_adapter_for_pe",
        "--failed_out",
        str(failed_reads),
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--low_complexity_filter",
        "--json",
        str(json),
        "--html",
        str(html),
    ]

    if len(input_reads) == 2 and len(output_reads) == 2:
        args += ["-I", str(input_reads[1]), "-O", str(output_reads[1])]
    elif len(input_reads) == 1 and len(output_reads) == 1:
        pass
    else:
        raise ValueError("Input and output reads must both be single or paired-end.")

    return sp.run(args, capture_output=True, text=True)


log_f = snakemake.log[0]  # type: ignore
with open(log_f, "w") as log:
    log.write("Initialized log and starting script...\n")
    result = None
    try:
        result = f(log)
    except BaseException as e:
        log.write(f"Error during run: {e}\n")
        log.write(traceback.format_exc())
        raise
    else:
        if result is not None:
            log.write(f"STDOUT: {result.stdout}\n")
            log.write(f"STDERR: {result.stderr}\n")
        log.write("Completed successfully.\n")
