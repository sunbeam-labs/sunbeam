import logging
from pathlib import Path


class ConditionalLevelFormatter(logging.Formatter):
    def format(self, record):
        # For WARNING and above, include "LEVELNAME: message"
        if record.levelno >= logging.WARNING:
            return f"{record.levelname}: {record.getMessage()}"
        # For lower levels like INFO, just return the message
        return record.getMessage()


def get_sunbeam_logger() -> logging.Logger:
    """Basic logger for general library output."""
    logger = logging.getLogger("sunbeam")
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(ConditionalLevelFormatter())
        logger.addHandler(ch)

    return logger


def get_pipeline_logger(log_file: Path = None) -> logging.Logger:
    """Sets up logging for the main pipeline entry point."""
    logger = logging.getLogger("sunbeam.pipeline")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if not logger.handlers:
        # Console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(ConditionalLevelFormatter())

        # File handler
        if log_file is None:
            raise ValueError(
                "log_file is None but the logger hasn't been initialized with a file handler yet"
            )
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))

        logger.addHandler(ch)
        logger.addHandler(fh)

    return logger


class ExtensionLoggerAdapter(logging.LoggerAdapter):
    def process(self, msg, kwargs):
        return f"[{self.extra['ext'].upper()}] {msg}", kwargs


def get_extension_logger(name: str) -> logging.Logger:
    """Returns a logger for a specific extension under the pipeline."""
    logger = logging.getLogger(f"sunbeam.pipeline.extensions.{name}")
    logger.setLevel(logging.DEBUG)

    # Remove any direct handlers — let it propagate to sunbeam.pipeline
    logger.handlers.clear()
    logger.propagate = True

    return ExtensionLoggerAdapter(logger, {"ext": name})
