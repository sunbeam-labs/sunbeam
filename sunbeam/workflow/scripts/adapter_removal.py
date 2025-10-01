import sys
import traceback

log_f = snakemake.log[0]  # type: ignore

with open(log_f, "w") as log:
    log.write("Starting adapter removal step...\n")
    result = None
    try:
        import subprocess as sp
        from pathlib import Path

        input_reads = snakemake.input.reads  # type: ignore
        output_reads = snakemake.output.reads  # type: ignore
        failed_reads = snakemake.output.failed  # type: ignore
        json = snakemake.output.json  # type: ignore
        adapter = snakemake.params.adapter  # type: ignore
        threads = snakemake.threads  # type: ignore

        adapter_fp = Path(adapter)
        assert adapter_fp.is_file(), f"Adapter file {adapter} does not exist."
        assert (
            adapter_fp.stat().st_size > 0
        ), f"Adapter file {adapter} is empty (use a dummy file like `> A\nA` to skip this step)."

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
            str(failed_reads[0]),
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--low_complexity_filter",
            "--json",
            str(json),
        ]

        if len(input_reads) == 2 and len(output_reads) == 2:
            args += ["-I", str(input_reads[1]), "-O", str(output_reads[1])]
        elif len(input_reads) == 1 and len(output_reads) == 1:
            pass
        else:
            raise ValueError(
                "Input and output reads must both be single or paired-end."
            )

        result = sp.run(args, capture_output=True, text=True)

        log.write("Adapter removal completed successfully.\n")
    except BaseException as e:
        log.write(f"Error during adapter removal: {e}\n")
        log.write(traceback.format_exc())

        if result is not None:
            captured = result.stderr
            if captured:
                log.write("Captured stderr before failure:\n")
                log.write(captured)

        raise
