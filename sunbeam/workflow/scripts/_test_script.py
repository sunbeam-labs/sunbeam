from pathlib import Path


l = snakemake.log[0]
o = Path(snakemake.output[0])

with open(l, "w") as log:
    log.write("HERE\n")

    from sunbeam.logging import get_pipeline_logger, get_extension_logger
    from sunbeam.project import SampleList, SunbeamConfig, output_subdir

    logger = get_pipeline_logger()
    logger.info("This works!")

    SampleList()

    o.write_text("This is a test output file.\n")
