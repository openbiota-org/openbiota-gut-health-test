"""Allow ``python -m openbiota``."""

from __future__ import annotations

import sys

from openbiota.cli import main

if __name__ == "__main__":
    sys.exit(main())
