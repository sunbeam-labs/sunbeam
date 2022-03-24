"""
Summarize stats for coverage data for each sample 
"""

import csv
import numpy

input_fp = snakemake.input[0]
sample = snakemake.wilcards.sample
output_fp = snakemake.output[0]

with open(input_fp) as f:
    reader = csv.reader(f, delimiter='\t')    
    data = {}
    for row in reader:
        if not data.get(row[0]):
            data[row[0]] = []
        data[row[0]].append(int(row[2]))

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

# write out stats per segment per sample
fields = ['sample','contig', 'sum', 'min', 'max', 'mean', 'median', 'stddev', 'coverage', 'length']
with open(output_fp, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
    writer.writerows(output_rows)