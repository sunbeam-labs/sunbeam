from typing import Dict, Iterator, Tuple, Union
from sunbeam.bfx.parse import parse_sam


def get_mapped_reads(fp: str, min_pct_id: float, min_len_frac: float) -> Iterator[str]:
    """
    Takes a SAM file and returns an iterator of read names that are mapped
    """
    with open(fp, "r") as sam_file:
        for read in parse_sam(sam_file):
            if (
                (not read["FLAG"] & 0x4)  # not unmapped
                and (_get_frac(read) > min_len_frac)
                and (_get_pct_identity(read) > min_pct_id)
            ):
                yield read["QNAME"]


def _get_pct_identity(
    read: Dict[str, Union[int, float, str, Tuple[int, str]]]
) -> float:
    edit_dist = read.get("NM", 0)
    pct_mm = float(edit_dist) / len(read["SEQ"])
    return 1 - pct_mm


def _get_frac(read: Dict[str, Union[int, float, str, Tuple[int, str]]]) -> float:
    clip = 0
    for pair in read["CIGAR"]:
        if pair[1] == "S" or pair[1] == "H":
            clip += pair[0]
    frac = float(len(read["SEQ"])) / (len(read["SEQ"]) + clip)
    return frac
