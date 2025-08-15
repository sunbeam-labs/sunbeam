import logging
from pathlib import Path


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
)


def get_sunbeam_logger() -> logging.Logger:
    return logging.getLogger()


def get_pipeline_logger(log_fp: Path | None = None) -> logging.Logger:
    return logging.getLogger()


def get_extension_logger(name: str) -> logging.Logger:
    return logging.getLogger()
