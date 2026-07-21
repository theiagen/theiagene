"""Enable ``python -m theiagene`` to invoke the unified CLI."""

import sys

from theiagene.cli import main


if __name__ == "__main__":
    sys.exit(main())
