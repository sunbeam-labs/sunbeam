from snakemake.utils import listfiles

def validate_paths(cfg):
    """Process paths in config file.

    Any keys ending in _fp are converted to Paths, any ~ are expanded,
    and all paths except for `output_fp` are checked to make sure they
    exist.
    
    :param: cfg a config file
    :param: required_paths a list of required paths
    :returns: an updated copy of cfg
    """
    new_cfg = dict()
    for k, v in cfg.items():
        if k.endswith('_fp'):
            v = Path(v).expanduser()
            if k != 'output_fp' and not v.exists():
                raise ValueError(
                    "Path for `%s` does not exist: %s" % (k,v))
        new_cfg[k] = v
    return new_cfg


def build_sample_list(data_fp, filename_fmt):
    """Build a list of samples from a data filepath and filename format.

    :param: data_fp a Path to data files
    :param: filename_fmt a string giving wildcards for {sample} and (optionally) 
    {rp}, for read pair (e.g. R1 or R2).
    :returns: A dictionary of samples, with sample names as keys
    """
    files = list(listfiles(str(data_fp/filename_fmt)))
    Samples = {t[1]['sample']: {} for t in files}
    for f in files:
        fpath = f[0]
        wcards = f[1]
        rp = wcards['rp'] if 'rp' in wcards.keys() else False
        if rp:
            Samples[wcards['sample']][rp] = fpath
            Samples[wcards['sample']]['paired'] = True
        else:
            Samples[wcards['sample']]['file'] = fpath
            Samples[wcards['sample']]['paired'] = False
    return Samples


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
