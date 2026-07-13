"""FunVIP exception hierarchy.

Every error FunVIP raises on purpose derives from :class:`FunVIPError`, so the
top-level handler in ``main`` (and any caller) can tell an *expected, explained*
FunVIP failure apart from an unexpected bug and report it cleanly.

Guidance:
- Raise the most specific subclass that fits, with a message that names the
  offending thing (group, gene, FI id, file path, option) -- never a bare
  ``raise Exception``.
- Reserve a plain ``FunVIPError`` for cases that do not fit a subclass.
- Do not use these for internal control flow; they mean "stop and tell the user".
"""


class FunVIPError(Exception):
    """Base class for all FunVIP errors."""


class ConfigError(FunVIPError):
    """Invalid or inconsistent command-line options / configuration."""


class InputError(FunVIPError):
    """Malformed or missing query/database input (sequences, tables, metadata)."""


class SearchError(FunVIPError):
    """Failure in the BLAST/mmseqs search stage or its result matrices."""


class ClusterError(FunVIPError):
    """Failure assigning sequences to groups or selecting an outgroup."""


class DatasetError(FunVIPError):
    """Failure building a group/gene dataset or its (concatenated) alignment."""


class TreeError(FunVIPError):
    """Failure building, rooting, or interpreting a phylogenetic tree."""


class ExternalToolError(FunVIPError):
    """An external tool (mafft, trimal, fasttree, iqtree, raxml, mmseqs, gblocks,
    modeltest, ...) was missing, failed, or returned an unusable result."""
