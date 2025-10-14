import traceback
from typing import TextIO


def f(log: TextIO):
    import pandas
    from sunbeam.bfx.reports import parse_fastqc_quality

    input_reports = snakemake.input.reports  # type: ignore
    output_report = snakemake.output.report  # type: ignore

    quality_list = [parse_fastqc_quality(file) for file in input_reports]
    quality_list = [qr for qr in quality_list if qr is not None]
    quality_table = pandas.concat(quality_list, axis=1).transpose()
    quality_table.to_csv(output_report, sep="\t", index_label="Samples")


log_f = snakemake.log[0]  # type: ignore
with open(log_f, "w") as log:
    log.write("Initialized log and starting script...\n")
    try:
        f(log)
    except BaseException as e:
        log.write(f"Error during run: {e}\n")
        log.write(traceback.format_exc())
        raise
    else:
        log.write("Completed successfully.\n")
