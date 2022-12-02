import csv
import datetime
import os

def parse_err_and_warn(log_fp: str) -> tuple:
    """Go through output log and return a list of warnings and a list of errors"""
    warns = list()
    errs = list()


    return warns, errs

def compile_benchmarks(benchmark_fp: str, stats_fp: str):
    """Aggregate all the benchmark files into one and put it in stats_fp"""
    benchmarks = os.listdir(benchmark_fp)
    if not benchmarks:
        print("No benchmark files found")
        return None
    headers = ['rule']
    with open(os.path.join(benchmark_fp, benchmarks[0]), 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers += next(reader)
    
    if not os.path.exists(stats_fp):
        os.makedirs(stats_fp)
    stats_file = os.path.join(stats_fp, f"{str(int(datetime.datetime.now().timestamp() * 1000))}_benchmarks.tsv")
    with open(stats_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        for fp in benchmarks:
            with open(os.path.join(benchmark_fp, fp), 'r') as g:
                reader = csv.reader(g, delimiter='\t')
                next(reader) # Headers line
                writer.writerow(fp[:-4] + next(reader))


