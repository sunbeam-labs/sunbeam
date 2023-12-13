"""
Supporting functions for QC rules.
"""

import gzip
from more_itertools import grouper
from pathlib import Path
from sunbeamlib.parse import write_many_fastq
from typing import List, TextIO


def filter_ids(fp_in: Path, fp_out: Path, ids: List[str], log: TextIO) -> None:
    """Remove ids from FASTQ file.
    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        if not ids:
            log.write("Komplexity IDs list empty\n")
        else:
            ids = set(ids)
            records_dict = dict()
            records = []
            print("making dictionary and records")
            for g in grouper(f_in.readlines(), 4):
                header_str = g[0][1:].strip()
                seq_str = g[1].strip()
                plus_str = g[2].strip()
                quality_str = g[3].strip()
                records.append((header_str, seq_str, plus_str, quality_str))
                newkey = header_str.split(" ")[0]
                if newkey[-2:] == "/1" or newkey[-2:] == "/2":
                    newkey = newkey[:-2]
                records_dict[newkey] = (header_str, seq_str, plus_str, quality_str)
            print("length of unfiltered dict")
            print(len(records_dict))
            print("length of unfiltered list")
            print(len(records))
            print("make set of ids")
            records_ids = [*records_dict]  # grab keys
            records_ids = set(records_ids)
            print(len(records_ids))
            print("length of unwanted ids")
            print(len(ids))
            log.write(
                f"Unfiltered read count: {len(records_ids)}\nKomplexity IDs to be filtered {len(ids)}\n"
            )
            if records_ids.issuperset(ids):  # records contains every element of ids
                filtered_records_ids = records_ids.difference(ids)
                print("Length of expected files")
                print(len(filtered_records_ids))
                #                from joblib import Parallel, delayed  # attempt to parallelize TODO
                #                Parallel(n_jobs=2)(delayed(records.remove(records_dict[ids[i]]))(i ** 2) for i in range(len(ids)))
                for unwanted_key in ids:
                    records.remove(records_dict[unwanted_key])
                print("Length of filtered file")
                print(len(records))
            elif records_ids.isdisjoint(ids):  # there are no shared items between sets
                log.write(
                    "None of the low complexity read IDs found in unfiltered reads file! Please fix me.\n"
                )
                raise ValueError(
                    "IDs provided should not be missing from the unfiltered fastq file"
                )


def filter_ids_stash(fp_in: Path, fp_out: Path, ids: List[str], log: TextIO) -> None:
    """Remove ids from FASTQ file.
    fp_in: path to input FASTQ
    fp_out: path to output FASTQ
    ids: list of ids to be removed
    """
    with gzip.open(fp_in, "rt") as f_in, gzip.open(fp_out, "wt") as f_out:
        records_dict = dict()
        records = []
        for g in grouper(f.readlines(), 4):
            header_str = g[0][1:].strip()
            seq_str = g[1].strip()
            plus_str = g[2].strip()
            quality_str = g[3].strip()
            records.extend((header_str, seq_str, plus_str, quality_str))
            newkey = header_str.split(" ")[0]
            records_dict[newkey] = (header_str, seq_str, plus_str, quality_str)
        print(len(records_dict))
        print(len(records))
        # for r in parse_fastq(f_in):
        #    newkey = r[0].split(' ')[0]
        #    records_dict[newkey] = r
        records_ids = [*records_dict]  # grab keys
        records_ids = set(records_ids)

        # try different parsing
        # read_ids = []
        # for g in grouper(f.readlines(), 4):
        #    read_ids.append(g[0][1:].strip().split(' ')[0])
        # read_ids = set(read_ids)
        # unfiltered = [r for r in parse_fastq(f_in)]

        ids = set(ids)
        log.write(
            f"Unfiltered read count: {len(records_ids)}\nKomplexity IDs to be filtered {len(ids)}\n"
        )
        if records_ids.issuperset(ids):  # records contains every element of ids
            filtered_records_ids = records_ids.difference(ids)
            for unwanted_key in filtered_records_ids:
                #                del records_dict[unwanted_key]
                records.remove(records_dict[unwanted_key])
            log.write(f"Filtered read count: {len(records_dict)}\n")
            #            write_many_fastq(list(records_dict.values()), f_out)
            write_many_fastq(records, f_out)
        elif records_ids.isdisjoint(ids):  # there are no shared items between sets
            log.write(
                "None of the low complexity read IDs found in unfiltered reads file! Please fix me.\n"
            )
            raise ValueError(
                "IDs provided should not be missing from the unfiltered fastq file"
            )


def remove_pair_id(id: str, log: TextIO) -> str:
    """Remove the 1 or 2 from a paired read ID

    id: id string
    """
    id = id.strip()
    if id[-2:] == "/1" or id[-2:] == "/2":
        return id[:-2]

    # Assuming it's the newer id variant where komplexity removes the second half (containing pair number)
    return id
