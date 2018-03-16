import sys
import ruamel.yaml
import argparse

from sunbeamlib import config

def main(argv):
    # Two subcommands: update and modify config files
    #   - update: updates an old config file to a new version (possibly in-place)
    #   - modify: alters keys in an existing config file with new values

    def add_common_args(parser):
        output_args = parser.add_mutually_exclusive_group()
        output_args.add_argument(
            "-i", "--in_place", action="store_true",
            help="Alters config file in place")
        output_args.add_argument(
            "-o", "--out", type=argparse.FileType('w'), metavar="FILE",
            help="Where to write modified config file", default=sys.stdout)
        parser.add_argument(
            "config_file", type=argparse.FileType('r'),
            help="Existing config file to update/modify")
        return(parser)
    
    parser = argparse.ArgumentParser("sunbeam config")
    parser.set_defaults(func=lambda x: parser.print_help())
    subcommands = parser.add_subparsers()

    # 'update' subcommand
    update_desc_str = (
        "Usage examples:\n"
        "1. To update a config file in place:\n"
        "    $ sunbeam config update -i my_config.yml\n"
        "2. To write an update copy to a new file:\n"
        "    $ sunbeam config update old_config.yml -o new_config.yml")
    update_command = subcommands.add_parser(
        'update', description=update_desc_str,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    update_command.set_defaults(func=update)
    update_command.add_argument(
        "-t", "--template", default=None, type=argparse.FileType('r'),
        metavar="FILE", help="Path to custom config file template, in YAML format")
    update_command.add_argument(
        "--strict", action="store_true",
        help="Remove keys that no longer exist in the new config file")
    update_command = add_common_args(update_command)    

    # 'modify' subcommand
    modify_desc_str = (
        "Usage examples:\n"
        "1. To apply a set of defaults to an existing config file in place:\n"
        "    $ sunbeam config modify -i -f defaults.yml my_config.yml\n"
        "2. To change a single key:value pair in the 'mapping' section:\n"
        "    $ sunbeam config modify -i -s 'mapping: {keep_unaligned: True}'")
    modify_command = subcommands.add_parser(
        "modify", description=modify_desc_str,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    modify_command.set_defaults(func=modify)
    modifier_args = modify_command.add_mutually_exclusive_group()
    modifier_args.add_argument(
        "-s", "--str", type=str, metavar="STR",
        help="YAML string (e.g. 'blast: {threads: 4}')")
    modifier_args.add_argument(
        "-f", "--file", type=argparse.FileType("r"), metavar="FILE",
        help="YAML file with new config values") 
    modify_command = add_common_args(modify_command)
    
    args = parser.parse_args(argv)
    args.func(args)
    
def update(args):
    if args.in_place and args.strict:
        raise SystemExit(
            "Incompatible arguments: --strict cannot be used with --in_place "
            "to prevent unexpected loss of config settings.")

    old_config = ruamel.yaml.safe_load(args.config_file)

    # Remove the old version number
    old_config.get('all', {}).pop('version', None)

    new_config = config.new(
        conda_fp = config._find_conda_fp(),
        project_fp = old_config.get('all', {}).get('version', ''),
        template = args.template)

    updated_config = config.update(
        new_config, old_config, args.strict)

    # If strict, preserve the blast databases
    if args.strict:
        blastdbs = old_config.get('blastdbs', None)
        updated_config['blastdbs'] = blastdbs
    
    if args.in_place:
        cfg_fp = args.config_file.name
        args.config_file.close()
        with open(cfg_fp, 'w') as out:
            config.dump(updated_config, out)
    else:
        config.dump(updated_config, args.out)

def modify(args):
    update_src = args.str if args.str else args.file
    new_values = ruamel.yaml.safe_load(update_src)
    if isinstance(new_values, str) and args.str:
        raise SystemExit(
            "Invalid YAML in --str. Did you make sure to put spaces between "
            "keys and values?")
    modified_config = config.update(args.config_file, new_values, strict=False)
    if args.in_place:
        cfg_fp = args.config_file.name
        args.config_file.close()
        with open(cfg_fp, 'w') as out:
            config.dump(modified_config, out)
    else:
        config.dump(modified_config, args.out)

    
        
