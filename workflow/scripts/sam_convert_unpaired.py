with open(snakemake.log[0], "w") as log:
    with open(snakemake.input[0]) as sam, open(snakemake.output[0], "w") as fq:
        for line in sam.readlines():
            arr = line.strip().split("\t")
            if line[0] == "@":
                continue
            fq.write(f"{arr[0]}\n{arr[9]}\n+\n{arr[10]}\n")
            