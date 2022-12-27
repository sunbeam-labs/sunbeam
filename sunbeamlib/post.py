import csv
import datetime
import os


WARN_STR = ['warn', 'warning']
ERR_STR = ['err', 'error', 'exception']


def parse_err_and_warn(log_fp: str) -> tuple:
    """Go through output log and return a list of warnings and a list of errors"""
    warns = list()
    errs = list()

    with open(log_fp) as f:
        for n, l in enumerate(f, 1):
            if [1 for w in WARN_STR if w in l.lower()]:
                warns.append(n)
            if [1 for e in ERR_STR if e in l.lower()]:
                errs.append(n)

    return warns, errs


def parse_rule_logs(log_fp: str) -> list:
    """Go through each per-rule log and report any warnings or errors"""
    log_fps = os.listdir(log_fp)
    alerts = []
    for fp in log_fps:
        warns, errs = parse_err_and_warn(os.path.join(log_fp, fp))
        if warns or errs:
            alerts.append(f"{fp}:\nWarnings: ({len(warns)}) {warns}\nErrors: ({len(errs)}) {errs}\n")
    
    return alerts


def compile_benchmarks(benchmark_fp: str, stats_fp: str):
    """Aggregate all the benchmark files into one and put it in stats_fp"""
    benchmarks = []
    try:
        benchmarks = os.listdir(benchmark_fp)
    except FileNotFoundError as e:
        print("No benchmark files found")
        return None
    if not benchmarks:
        print("No benchmark files found")
        return None
    headers = ["rule"]
    with open(os.path.join(benchmark_fp, benchmarks[0]), "r") as f:
        reader = csv.reader(f, delimiter="\t")
        headers += next(reader)

    if not os.path.exists(stats_fp):
        os.makedirs(stats_fp)
    stats_file = os.path.join(
        stats_fp,
        f"{str(int(datetime.datetime.now().timestamp() * 1000))}_benchmarks.tsv",
    )
    with open(stats_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for fp in sorted(benchmarks):
            with open(os.path.join(benchmark_fp, fp), "r") as g:
                reader = csv.reader(g, delimiter="\t")
                next(reader)  # Headers line
                writer.writerow([fp[:-4]] + next(reader))
