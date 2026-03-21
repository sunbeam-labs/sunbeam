import pathlib
import csv
import json
import math
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


def rows_from_dicts(ds, empty_value=""):
    # assume that first dict has all keys
    # use order of keys in first dict
    keys = None
    for d in ds:
        if keys is None:
            keys = list(d.keys())
            yield keys
        vals = [d.get(k, empty_value) for k in keys]
        yield vals


def make_preprocess_report(
    adapter_fps, trim_fps, complexity_fps, decontam_fps, output_fp, samples, log_fp
):
    per_sample_fps = zip(samples, adapter_fps, trim_fps, complexity_fps, decontam_fps)
    per_sample_data = []
    for sample, adapt_fp, trim_fp, comp_fp, deco_fp in per_sample_fps:
        row = {"sample": sample}

        adapter_data = parse_adapter_report(adapt_fp)
        row.update(adapter_data)
        trim_data = parse_quality_report(trim_fp)
        row.update(trim_data)
        complexity_data = parse_complexity_report(comp_fp)
        row.update(complexity_data)
        decontam_data = parse_decontam_report(deco_fp)
        row.update(decontam_data)

        per_sample_data.append(row)

    rows = rows_from_dicts(per_sample_data)
    with open(output_fp, "w") as f:
        for row in rows:
            f.write("\t".join(str(val) for val in row))
            f.write("\n")


def make_fastqc_report(
    input_report_fps,
    output_report_fp,
    log_fp,
):
    qual_data_long = []
    for fp in input_report_fps:
        fp = pathlib.Path(fp)
        sample_id = (fp.resolve().parent.name).split("_fastqc")[0]  # Hate this
        with open(fp) as f:
            for base_idx, mean_qual in parse_fastqc_quality(f):
                qual_data_long.append((sample_id, base_idx, mean_qual))
    with open(output_report_fp, "w") as f:
        for row in simple_pivot(qual_data_long, "Samples", ""):
            f.write("\t".join(row))
            f.write("\n")


def simple_pivot(triples, id_colname, empty_val):
    # We need to iterate through triples twice
    triples = list(triples)

    colnames = {}
    colnames[id_colname] = None
    for _, colname, _ in triples:
        colnames[colname] = None
    yield list(colnames)

    nrow = 0
    current_rowname = None
    current_vals = {}
    for rowname, colname, val in triples:
        if rowname != current_rowname:
            # Emit previous row
            if nrow > 0:
                yield list(current_vals.get(k, empty_val) for k in colnames)
            # Set up current row
            current_rowname = rowname
            current_vals = {}
            current_vals[id_colname] = rowname

        current_vals[colname] = val
        nrow += 1

    # Emit final row if there were any
    if nrow > 0:
        yield list(current_vals.get(k, empty_val) for k in colnames)


def parse_fastqc_quality(f):
    in_per_base_quals = False
    header_seen = False
    for line in f:
        if in_per_base_quals:
            line = line.rstrip()
            if not header_seen:
                assert line.startswith("#Base")
                header_seen = True
            elif line.startswith(">>END_MODULE"):
                break
            else:
                vals = line.split("\t")
                yield (vals[0], vals[1])
        elif line.startswith(">>Per base sequence quality"):
            in_per_base_quals = True
