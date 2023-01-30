import pysam
import sys


def get_mapped_reads(fp, min_pct_id, min_len_frac, log):
    log.write(
        f"Running get_mapped_reads on {fp} with min_pct_id = {min_pct_id} and min_len_frac = {min_len_frac}\n"
    )
    sam = pysam.AlignmentFile(fp)
    count = 0
    count_ret = 0
    for read in sam:
        count += 1
        log.write(str(read.is_unmapped))
        if (
            (not read.is_unmapped)
            and (_get_frac(read) > min_len_frac)
            and (_get_pct_identity(read) > min_pct_id)
        ):
            count_ret += 1
            yield read.query_name

    log.write(f"Total reads processed: {count}")
    log.write(f"Total reads used: {count_ret}")


def _get_pct_identity(read):
    if read.has_tag("NM"):
        edit_dist = read.get_tag("NM")
    else:
        edit_dist = 0
    pct_mm = float(edit_dist) / read.alen
    return 1 - pct_mm


def _get_frac(read):
    cigar = read.cigartuples
    clip = 0
    for pair in cigar:
        if pair[0] == 4 or pair[0] == 5:
            clip = clip + pair[1]
    frac = float(read.query_alignment_length) / (read.query_alignment_length + clip)
    return frac
