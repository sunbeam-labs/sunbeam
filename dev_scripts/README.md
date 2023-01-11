# dev_scripts

The utilities in this directory are aimed at sunbeam developers. Run these scripts from the main sunbeam directory, e.g.

```
bash dev_scripts/reformat.sh
```

The sunbeam_config.yml file is included to run some snakemake commands.

### generate_archive.sh

Generate a tarball for a sunbeam release. This requires first running snakemake's `--archive` to generate the `.snakemake/conda-archive` directory and then manually wrapping up the necessary directories for release.

### generate_dockerfile.sh

This is for using singularity - a feature that still isn't implemented in sunbeam. If this is implemented it will have a similar release structure as archives.

### reformat.sh

This is used for formatting all python and snakemake files using `black` and `snakefmt`.