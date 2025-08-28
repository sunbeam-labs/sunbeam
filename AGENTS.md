# Dear Agents,

This is a Python library wrapping a Snakemake pipeline for metagenomic bioinformatics.

## Project Structure

- `sunbeam/`
  - `bfx/`: Helper bioinformatics functions
  - `configs/`: Example and template config and profile files
  - `project/`: Classes for managing Sunbeam specific config files
  - `scripts/`: CLI commands
  - `workflow/`: The Snakemake workflow definition
- `extensions/`: This is where extensions are installed to (`git clone`d). This directory can also live elsewhere as defined by the `$SUNBEAM_EXTENSIONS` variable

## Contribution Guidelines

Before making contributions to the codebase:

- Reformat the codebase with `black .`
- Make sure unit tests pass with `pytest tests/unit/`
  - If they aren't passing be sure to include an explanation why
