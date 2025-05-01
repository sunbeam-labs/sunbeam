import gzip
import os
import shutil
import subprocess as sp
import sys


with open(snakemake.log[0], "w") as log:
    fwd_adapters = snakemake.config["qc"]["fwd_adapters"]
    rev_adapters = snakemake.config["qc"]["rev_adapters"]
    if fwd_adapters or rev_adapters:
        overlap = float("inf")
        if fwd_adapters:
            overlap = min(min(len(a) for a in fwd_adapters), overlap)
            fwd_adapter_str = "-b " + " -b ".join(
                snakemake.config["qc"]["fwd_adapters"]
            )
        if rev_adapters:
            overlap = min(min(len(a) for a in rev_adapters), overlap)
            rev_adapter_str = "-B " + " -B ".join(
                snakemake.config["qc"]["rev_adapters"]
            )

        try:
            args = [
                "cutadapt",
                "-O",
                str(overlap),
                "--cores",
                str(snakemake.threads),
            ]
            args += (
                snakemake.config["qc"]["cutadapt_opts"].split(" ")
                if snakemake.config["qc"]["cutadapt_opts"]
                else []
            )
            args += fwd_adapter_str.split(" ")
            args += rev_adapter_str.split(" ")
            args += [
                "-o",
                f"{snakemake.output.ngz1}",
                "-p",
                f"{snakemake.output.ngz2}",
                f"{snakemake.input.r1}",
                f"{snakemake.input.r2}",
            ]
            cutadapt_output = sp.check_output(
                args,
                stderr=sp.STDOUT,
            )
        except sp.CalledProcessError as e:
            log.write(e.output.decode())
            sys.exit(e.returncode)
        log.write(cutadapt_output.decode())

        with (
            open(snakemake.output.ngz1, "rb") as f_in,
            gzip.open(snakemake.output.r1, "wb") as f_out,
        ):
            shutil.copyfileobj(f_in, f_out)
        with (
            open(snakemake.output.ngz2, "rb") as f_in,
            gzip.open(snakemake.output.r2, "wb") as f_out,
        ):
            shutil.copyfileobj(f_in, f_out)
    else:
        log.write("Adapters not found, skipping adapter removal...")
        os.symlink(snakemake.input.r1, snakemake.output.r1)
        os.symlink(snakemake.input.r2, snakemake.output.r2)
