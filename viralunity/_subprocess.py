"""Shared subprocess helper for ViralUnity CLI commands.

Internal module — not part of the public API. Used by the `get-databases`
and `build-deacon-index` subcommands to invoke external tools (`wget`,
`tar`, `diamond`, `deacon`, `datasets`, ...) with consistent echoing and
error handling.
"""

from __future__ import annotations

import subprocess
from typing import Optional, Sequence

import click


def run_command(cmd: Sequence[object], cwd: Optional[str] = None) -> None:
    """Run a subprocess command, streaming output and raising on failure.

    Args:
        cmd: Argument list (no shell). Each element is coerced to ``str``
            when echoed; pass paths and flags as separate elements.
        cwd: Optional working directory in which to run the command.

    Raises:
        click.ClickException: If the subprocess exits non-zero. The exit
            code and command are included in the message so the failure
            surfaces cleanly in the CLI output.
    """
    click.echo(f"$ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(list(cmd), cwd=cwd)
    if result.returncode != 0:
        raise click.ClickException(
            f"Command failed with exit code {result.returncode}: "
            f"{' '.join(str(c) for c in cmd)}"
        )
