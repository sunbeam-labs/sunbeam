def index_files(genome):
    """Return the bowtie index files for a file."""
    fwd, rev = (
        expand(
            "{index_fp}/{genome}.{index}.bt2",
            index_fp=Cfg['bt2_index_fp'],
            genome=genome,
            index=range(1,5)),
        expand(
            "{index_fp}/{genome}.rev.{index}.bt2",
            index_fp=Cfg['bt2_index_fp'],
            genome=genome,
            index=range(1,3))
    )
    return fwd + rev
