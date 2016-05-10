from pathlib import Path

import yaml

from snakemake.workflow import expand
from snakemake.utils import listfiles

def makepath(path):
    return Path(path).expanduser()


def verify(path):
    if path.exists():
        return path.resolve()
    else:
        raise ValueError("Path %s does not exist" % path)

def validate_paths(cfg):
    """Process paths in config file.

    Any keys ending in _fp are converted to Paths, any ~ are expanded,
    and all paths except for `output_fp` are checked to make sure they
    exist.
    
    :param cfg: a config file
    :returns: an updated copy of cfg
    """
    new_cfg = dict()
    for k, v in cfg.items():
        if k.endswith('_fp'):
            v = makepath(v)
            print(v)
            if not v.is_absolute():
                raise ValueError("All paths must be absolute: %s:%s" %
                                 (k, str(v)))
            if k != 'output_fp':
                verify(v)
        new_cfg[k] = v
    return new_cfg

def load_subconfig(fp, verify=True):
    if verify:
        return validate_paths(yaml.load(open(fp)))
    else:
        return yaml.load(open(fp))

def process_databases(db_dict):
    """Process the list of databases.

    Expands the nucleotide and protein databases specified
    """
    dbs = {'nucl':{}, 'prot':{}}
    root = verify(makepath(db_dict['root_fp']))
    nucl = db_dict.get('nucleotide')
    prot = db_dict.get('protein')
    if nucl:
        dbs['nucl'] = {db: str(root/path) for db, path in nucl.items()}
    if prot:
        dbs['prot'] = {db: str(root/path) for db, path in prot.items()}
    return dbs
    
        
    
