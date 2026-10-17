"""Check VASP input files for common mistakes before you submit."""

__version__ = "0.1.0"

from .inputs import load
from .rules import RULES, Finding, run


def lint(directory, ignore=()):
    """Findings for one calculation directory."""
    return run(load(directory), ignore)


__all__ = ["RULES", "Finding", "__version__", "lint", "load", "run"]
