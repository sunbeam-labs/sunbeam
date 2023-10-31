import csv
import datetime
import os


def compile_benchmarks(benchmark_fp: str, stats_fp: str) -> None:
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

import sys
def update_benchmark(func):
    # added arguments inside the inner1,
    # if function takes any arguments,
    # can be added like this.
    def inner1(*args, **kwargs):
 
        # storing time before function execution
        sys.stderr.write("B")
         
        func(*args, **kwargs)
 
        # storing time after function execution
        sys.stderr.write("A")
 
    return inner1
