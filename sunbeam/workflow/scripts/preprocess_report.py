import traceback
from typing import TextIO


def f(log: TextIO):
    import pandas
    from sunbeam.bfx.reports import summarize_qual_decontam

    adapter_reports = snakemake.input.adapter  # type: ignore
    trim_reports = snakemake.input.trim  # type: ignore
    complexity_reports = snakemake.input.complexity  # type: ignore
    decontam_reports = snakemake.input.decontam  # type: ignore
    output_report = snakemake.output.report  # type: ignore
    samples = snakemake.params.samples  # type: ignore

    reports = {
        s: {"adapter": a, "trim": t, "complexity": c, "decontam": d}
        for s, a, t, c, d in zip(
            samples, adapter_reports, trim_reports, complexity_reports, decontam_reports
        )
    }
    log.write(f"Reports: {reports}\n")

    summary_list = [
        summarize_qual_decontam(q, d, k, paired_end)
        for q, d, k in zip(
            snakemake.input.trim_files,
            snakemake.input.decontam_files,
            snakemake.input.komplexity_files,
        )
    ]
    _reports = pandas.concat(summary_list)
    _reports.to_csv(snakemake.output[0], sep="\t", index_label="Samples")


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
