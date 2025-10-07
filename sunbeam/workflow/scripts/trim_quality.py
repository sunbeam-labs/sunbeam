import traceback
from contextlib import redirect_stderr
from io import StringIO
from typing import TextIO


def f(log: TextIO):
    from heyfastqlib.command import heyfastq_main

    input_reads = snakemake.input.reads  # type: ignore
    output_reads = snakemake.output.reads  # type: ignore
    output_count = snakemake.output.counter  # type: ignore
    window_width, window_threshold = snakemake.params.window  # type: ignore
    start_threshold = snakemake.params.start_threshold  # type: ignore
    end_threshold = snakemake.params.end_threshold  # type: ignore
    min_length = snakemake.params.min_length  # type: ignore

    stderr_capture = StringIO()

    # Redirect stderr to the buffer and call the function
    with redirect_stderr(stderr_capture):
        heyfastq_main(
            [
                "trim-qual",
                "--input",
                *input_reads,
                "--output",
                *output_reads,
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
            ]
        )

    # Retrieve the captured stderr output
    captured_stderr = stderr_capture.getvalue()

    with open(output_count, "w") as count:
        log.write(f"Counts: {captured_stderr}\n")
        count.write(captured_stderr)


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
