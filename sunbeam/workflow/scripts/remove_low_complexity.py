import sys
import traceback

log_f = snakemake.log[0]  # type: ignore

with open(log_f, "w") as log:
    log.write("Starting low-complexity read removal step...\n")
    stderr_capture = None
    try:
        from contextlib import redirect_stderr
        from io import StringIO
        from heyfastqlib.command import heyfastq_main

        input_reads = snakemake.input.reads  # type: ignore
        output_reads = snakemake.output.reads  # type: ignore
        output_count = snakemake.output.counter  # type: ignore
        kmer_size = snakemake.params.kmer_size  # type: ignore
        min_kscore = snakemake.params.min_kscore  # type: ignore

        stderr_capture = StringIO()

        # Redirect stderr to the buffer and call the function
        with redirect_stderr(stderr_capture):
            heyfastq_main(
                [
                    "filter-kscore",
                    "--input",
                    *input_reads,
                    "--output",
                    *output_reads,
                    "--kmer-size",
                    str(kmer_size),
                    "--min-kscore",
                    str(min_kscore),
                ]
            )

        # Retrieve the captured stderr output
        captured_stderr = stderr_capture.getvalue()

        with open(output_count, "w") as count:
            count.write(captured_stderr)

        log.write("Low-complexity read removal completed successfully.\n")
    except BaseException as e:
        log.write(f"Error during low-complexity read removal: {e}\n")
        log.write(traceback.format_exc())

        if stderr_capture is not None:
            captured = stderr_capture.getvalue()
            if captured:
                log.write("Captured stderr before failure:\n")
                log.write(captured)

        raise
