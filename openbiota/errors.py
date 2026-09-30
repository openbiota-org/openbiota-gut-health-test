"""Exception hierarchy.

Every failure mode gets a specific type so the CLI can print an actionable
message instead of a traceback. Per the build spec: fail loudly, never
silently produce a number from partial data.
"""

from __future__ import annotations


class OpenBiotaError(Exception):
    """Base class for all expected, user-facing failures."""


class DependencyError(OpenBiotaError):
    """An external tool (DIAMOND, ...) is missing or unusable."""


class PanelError(OpenBiotaError):
    """A panel YAML file is missing, malformed, or semantically invalid."""


class ReferenceFetchError(OpenBiotaError):
    """A reference sequence fetch failed or returned an unusable result."""


class SearchError(OpenBiotaError):
    """The DIAMOND search failed."""


class InputError(OpenBiotaError):
    """Input FASTQ files are missing, unreadable, or inconsistent."""
