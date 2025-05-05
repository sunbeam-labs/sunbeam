import csv
import datetime
import os
from pathlib import Path
from snakemake.common import Rules
from typing import Dict


def compile_benchmarks(
    benchmark_fp: str, Cfg: Dict[str, Dict | str], rules: Rules
) -> None:
    """Aggregate all the benchmark files into one and put it in stats_fp"""
    stats_fp = Path(Cfg["all"]["root"]) / "stats"
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

    dt = str(int(datetime.datetime.now().timestamp() * 1000))
    stats_file = os.path.join(stats_fp, f"{dt}_benchmarks.tsv")

    with open(stats_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)
        for fp in sorted(benchmarks):
            with open(os.path.join(benchmark_fp, fp), "r") as g:
                reader = csv.reader(g, delimiter="\t")
                next(reader)  # Headers line
                writer.writerow([fp[:-4]] + next(reader))

    compile_file_stats(stats_fp, Cfg, dt, rules)


def compile_file_stats(
    stats_fp: str, Cfg: Dict[str, Dict | str], dt: str, rules: Rules
) -> None:
    """Collect data on all inputs and outputs (as long as they still exist at this point) as well as dbs"""
    file_stats = {}
    output_fp = Path(Cfg["all"]["root"]) / Cfg["all"]["output_fp"]

    stats_file = os.path.join(stats_fp, f"{dt}_file_stats.tsv")

    # Collect info on all files in params
    param_fps = set()
    for n, r in rules._rules.items():
        if r.params:
            for p in r.params:
                if isinstance(p, dict):
                    for k, v in p.items():
                        try:
                            param_fps.add(Path(v))
                        except TypeError:
                            continue
                elif isinstance(p, list):
                    for i in p:
                        try:
                            param_fps.add(Path(i))
                        except TypeError:
                            continue
                else:
                    try:
                        param_fps.add(Path(p))
                    except TypeError:
                        continue

    for fp in param_fps:
        if fp.exists():
            if fp.is_file():
                file_stats[fp] = fp.stat().st_size
            elif fp.is_dir():
                for file in os.listdir(fp):
                    nfp = fp / file
                    if nfp.exists() and nfp.is_file():
                        file_stats[nfp] = nfp.stat().st_size

    # Collect info on all files in output_fp
    for root, dirs, files in os.walk(output_fp):
        for file in files:
            fp = Path(root) / file
            file_stats[fp] = fp.stat().st_size

    # Collect info on all dbs in Cfg
    for section, values in Cfg.items():
        if section == "all":
            continue
        for k, v in values.items():
            try:
                vp = Path(v)
            except TypeError:
                continue
            if vp.exists():
                if vp.is_file():
                    file_stats[vp] = vp.stat().st_size
                elif vp.is_dir():
                    for file in os.listdir(vp):
                        fp = vp / file
                        file_stats[fp] = fp.stat().st_size

    # Write to file
    with open(stats_file, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["file", "size"])
        for fp, size in file_stats.items():
            writer.writerow([fp, size])
