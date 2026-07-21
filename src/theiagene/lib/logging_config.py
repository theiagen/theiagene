"""Shared logging configuration for theiagene commands."""

import logging

LOG_FORMAT = "%(asctime)s - %(levelname)s - %(message)s"


def configure_logging(level: int = logging.DEBUG) -> None:
    """Configure root logging once for a command invocation.

    The subcommand modules obtain their own ``logging.getLogger(__name__)``
    loggers; this only installs a handler/format on the root logger."""
    logging.basicConfig(level=level, format=LOG_FORMAT)
