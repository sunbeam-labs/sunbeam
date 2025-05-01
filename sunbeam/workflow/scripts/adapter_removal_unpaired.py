import gzip
import os
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
            rev_adapter_str = "-g " + " -g ".join(
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
                snakemake.output.ngz,
                snakemake.input[0],
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
            open(snakemake.output.ngz) as f_in,
            gzip.open(snakemake.output.r, "wt") as f_out,
        ):
            f_out.writelines(f_in.readlines())
    else:
        log.write("Adapters not found, skipping adapter removal...")
        os.symlink(snakemake.input[0], snakemake.output.r)
