import warnings
import yaml
from pathlib import Path

from snakemake.workflow import expand
from snakemake.utils import listfiles

def makepath(path):
    return Path(path).expanduser()


def verify(path):
    path = Path(path)
    if path.exists():
        return path.resolve()
    else:
        raise ValueError("Path %s does not exist" % path)
    
    
def validate_paths(cfg, root):
    """Process paths in config file subsection.

    For each key ending in _fp, the value is:
    - converted to a pathlib.Path
    - ensured to be an absolute path, by appending `root` if needed
    - ensured to exist, if it is not the value from `output_fp`
    - expanded home directory ~
    
    :param cfg: a config file subsection
    :returns: an updated copy of cfg
    """
    new_cfg = dict()
    for k, v in cfg.items():
        if k.endswith('_fp'):
            v = makepath(v)
            if not v.is_absolute():
                v = root/v
            if k != 'output_fp':
                try: v = verify(v)
                except ValueError:
                    raise ValueError(
                        "For key '%s': path '%s' does not exist" % (k,v))
        new_cfg[k] = v
    return new_cfg


def check_config(cfg):
    """Resolve root in config file, then validate paths."""
    # Remove HOME_DIR from further path validation 
    cfg.pop("HOME_DIR", None)

    if 'root' in cfg['all']:
        root = verify(cfg['all']['root'])
    else:
        root = Path.cwd()
    # Iteratively check paths for each subsection
    new_cfg = dict()
    for section, values in cfg.items():
        new_cfg[section] = validate_paths(values, root)
    return new_cfg


def output_subdir(cfg, section):
    return cfg['all']['output_fp']/cfg[section]['suffix']


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
    
        
    
