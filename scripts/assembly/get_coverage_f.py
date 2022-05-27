import csv
import numpy
from io import TextIOWrapper

def parse_depth(f: TextIOWrapper) -> dict:
    reader = csv.reader(f, delimiter='\t')    
    data = {}
    for row in reader:
        if not data.get(row[0]):
            data[row[0]] = []
        data[row[0]].append(int(row[2]))
    return data

def get_cov_stats(data: dict, sample: str) -> list:
    # summarize stats for all segments present and append to output
    output_rows = []
    for segment in data.keys():
        sumval     = numpy.sum(data[segment])
        minval     = numpy.min(data[segment])
        maxval     = numpy.max(data[segment])
        mean       = numpy.mean(data[segment])
        median     = numpy.median(data[segment])
        stddev     = numpy.std(data[segment])
        gen_cov    = len(list(filter(lambda x: x!=0, data[segment])))
        gen_length = len(data[segment])
        row = [sample, segment, sumval, minval, maxval, mean, median, stddev, gen_cov, gen_length]
        output_rows.append(row)
    return output_rows

def write_csv(f: TextIOWrapper, cov: list):
    # write out stats per segment per sample
    fields = ['sample','contig', 'sum', 'min', 'max', 'mean', 'median', 'stddev', 'coverage', 'length']
    writer = csv.writer(f)
    writer.writerow(fields)
    writer.writerows(cov)