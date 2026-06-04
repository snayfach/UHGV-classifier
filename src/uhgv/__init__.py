from importlib.metadata import version

from uhgv.subcommands import classify, download_database

__all__ = ["classify", "download_database", "__version__"]

__version__ = version("uhgv")
