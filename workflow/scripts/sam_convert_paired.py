with open(snakemake.log[0], "w") as log:
    with open(snakemake.input[0]) as sam, open(snakemake.output.rp_1, "w") as rp1, open(snakemake.output.rp_2, "w") as rp2:
        for line in sam.readlines():
            arr = line.strip().split("\t")
            if line[0] == "@":
                continue
            if arr[1] == "77":
                rp1.write(f"{arr[0]}\n{arr[9]}\n+\n{arr[10]}\n")
            elif arr[1] == "141":
                rp2.write(f"{arr[0]}\n{arr[9]}\n+\n{arr[10]}\n")
            else:
                log.write(f"Skipping read {arr[0]} with flag {arr[1]}\n")
            