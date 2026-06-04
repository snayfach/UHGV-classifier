from importlib.metadata import version

from uhgv.modules import classify, download

__all__ = ["classify", "download", "__version__"]

__version__ = version("uhgv")
