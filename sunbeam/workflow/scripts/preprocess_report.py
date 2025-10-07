import csv
import json
import math
import traceback
from typing import (
    Dict,
    Iterable as TypingIterable,
    List,
    Mapping,
    MutableMapping,
    Optional,
    TextIO,
    Union,
)


Number = Union[int, float]


def _ensure_list(value: Union[str, TypingIterable[str]]) -> List[str]:
    """Return *value* as a list of paths."""

    if isinstance(value, (list, tuple, set)):
        return [str(v) for v in value]
    return [str(value)]


def parse_adapter_report(
    report: Union[str, TypingIterable[str]],
) -> Dict[str, Optional[Number]]:
    """Parse adapter removal (fastp) JSON reports.

    The fastp report already combines read 1 and read 2 for a sample. However,
    the helper gracefully supports being passed multiple reports by summing the
    count-based fields and computing weighted averages where appropriate.
    """

    totals: Dict[str, Optional[Number]] = {
        "adapter_total_reads_before": 0,
        "adapter_total_bases_before": 0,
        "adapter_total_reads_after": 0,
        "adapter_total_bases_after": 0,
        "adapter_passed_filter_reads": 0,
        "adapter_low_quality_reads": 0,
        "adapter_too_many_N_reads": 0,
        "adapter_low_complexity_reads": 0,
        "adapter_too_short_reads": 0,
        "adapter_too_long_reads": 0,
        "adapter_duplication_rate": 0.0,
        "adapter_read1_mean_length": 0.0,
        "adapter_read2_mean_length": 0.0,
        "adapter_gc_content_before": 0.0,
        "adapter_gc_content_after": 0.0,
    }

    total_before_bases = 0.0
    total_after_bases = 0.0

    reports = _ensure_list(report)
    if not reports:
        return totals

    for path in reports:
        with open(path) as fh:
            data = json.load(fh)

        summary = data.get("summary", {})
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})
        duplication = data.get("duplication", {})
        filtering_result = data.get("filtering_result", {})

        totals["adapter_total_reads_before"] += before.get("total_reads", 0)
        totals["adapter_total_bases_before"] += before.get("total_bases", 0)
        totals["adapter_total_reads_after"] += after.get("total_reads", 0)
        totals["adapter_total_bases_after"] += after.get("total_bases", 0)

        totals["adapter_passed_filter_reads"] += filtering_result.get(
            "passed_filter_reads", 0
        )
        totals["adapter_low_quality_reads"] += filtering_result.get(
            "low_quality_reads", 0
        )
        totals["adapter_too_many_N_reads"] += filtering_result.get(
            "too_many_N_reads", 0
        )
        totals["adapter_low_complexity_reads"] += filtering_result.get(
            "low_complexity_reads", 0
        )
        totals["adapter_too_short_reads"] += filtering_result.get("too_short_reads", 0)
        totals["adapter_too_long_reads"] += filtering_result.get("too_long_reads", 0)

        duplication_rate = duplication.get("rate")
        if duplication_rate is not None:
            totals["adapter_duplication_rate"] += duplication_rate

        for key, field in (
            ("adapter_read1_mean_length", "read1_mean_length"),
            ("adapter_read2_mean_length", "read2_mean_length"),
        ):
            value = after.get(field)
            if value is not None:
                totals[key] = (totals[key] or 0.0) + value

        before_gc = before.get("gc_content")
        if before_gc is not None:
            total_before_bases += before.get("total_bases", 0) * before_gc
        after_gc = after.get("gc_content")
        if after_gc is not None:
            total_after_bases += after.get("total_bases", 0) * after_gc

    num_reports = len(reports)
    if num_reports > 0:
        if totals["adapter_duplication_rate"] is not None:
            totals["adapter_duplication_rate"] /= num_reports

        for key in ("adapter_read1_mean_length", "adapter_read2_mean_length"):
            if totals[key] is not None:
                totals[key] /= num_reports

    before_bases = totals["adapter_total_bases_before"]
    totals["adapter_gc_content_before"] = (
        total_before_bases / before_bases if before_bases else 0.0
    )

    after_bases = totals["adapter_total_bases_after"]
    totals["adapter_gc_content_after"] = (
        total_after_bases / after_bases if after_bases else 0.0
    )

    return totals


def _parse_read_count_line(line: str) -> Optional[Dict[str, Number]]:
    line = line.strip()
    if not line:
        return None

    # Example: "Input reads: total=4000, average_length=70.00"
    if not line.startswith(("Input reads", "Output reads")):
        return None

    _, payload = line.split(":", 1)
    parts = payload.strip().split(",")
    values: Dict[str, Number] = {}
    for part in parts:
        if "=" not in part:
            continue
        key, value = part.strip().split("=", 1)
        if key == "total":
            values[key] = int(value)
        elif key == "average_length":
            values[key] = float(value)
    return values if values else None


def _aggregate_read_count_reports(
    report: Union[str, TypingIterable[str]], prefix: str
) -> Dict[str, Optional[Number]]:
    totals = {
        f"{prefix}_input_reads": 0,
        f"{prefix}_output_reads": 0,
        f"{prefix}_input_average_length": math.nan,
        f"{prefix}_output_average_length": math.nan,
    }

    input_bases = 0.0
    output_bases = 0.0

    for path in _ensure_list(report):
        with open(path) as fh:
            for line in fh:
                parsed = _parse_read_count_line(line)
                if not parsed:
                    continue
                if line.startswith("Input reads"):
                    totals[f"{prefix}_input_reads"] += int(parsed.get("total", 0))
                    input_bases += parsed.get("total", 0) * parsed.get(
                        "average_length", 0.0
                    )
                elif line.startswith("Output reads"):
                    totals[f"{prefix}_output_reads"] += int(parsed.get("total", 0))
                    output_bases += parsed.get("total", 0) * parsed.get(
                        "average_length", 0.0
                    )

    if totals[f"{prefix}_input_reads"]:
        totals[f"{prefix}_input_average_length"] = (
            input_bases / totals[f"{prefix}_input_reads"]
        )
    if totals[f"{prefix}_output_reads"]:
        totals[f"{prefix}_output_average_length"] = (
            output_bases / totals[f"{prefix}_output_reads"]
        )

    return totals


def parse_quality_report(
    report: Union[str, TypingIterable[str]],
) -> Dict[str, Optional[Number]]:
    return _aggregate_read_count_reports(report, "quality")


def parse_complexity_report(
    report: Union[str, TypingIterable[str]],
) -> Dict[str, Optional[Number]]:
    return _aggregate_read_count_reports(report, "complexity")


def parse_decontam_report(
    report: Union[str, TypingIterable[str]],
) -> Dict[str, Optional[Number]]:
    totals: MutableMapping[str, Number] = {}
    for path in _ensure_list(report):
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                for key, value in row.items():
                    try:
                        numeric_value: Number = int(value)
                    except (TypeError, ValueError):
                        try:
                            numeric_value = float(value)
                        except (TypeError, ValueError):
                            continue
                    totals[f"decontam_{key}"] = (
                        totals.get(f"decontam_{key}", 0) + numeric_value
                    )

    return dict(totals)


def parse_reports(reports: Mapping[str, Mapping[str, Union[str, TypingIterable[str]]]]):
    import pandas

    rows = []
    for sample, paths in sorted(reports.items()):
        row: Dict[str, Optional[Number]] = {"sample": sample}

        adapter_data = parse_adapter_report(paths["adapter"])
        quality_data = parse_quality_report(paths["trim"])
        complexity_data = parse_complexity_report(paths["complexity"])
        decontam_data = parse_decontam_report(paths["decontam"])

        row.update(adapter_data)
        row.update(quality_data)
        row.update(complexity_data)
        row.update(decontam_data)
        rows.append(row)

    df = pandas.DataFrame(rows).set_index("sample")
    return df.sort_index()


def f(log: TextIO):
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

    df = parse_reports(reports)
    df.to_csv(output_report, sep="\t", index_label="sample")


def run_from_snakemake():  # pragma: no cover - executed within snakemake
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


if "snakemake" in globals():  # pragma: no cover - executed within snakemake
    run_from_snakemake()
