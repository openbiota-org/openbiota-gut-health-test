#!/usr/bin/env python3
"""Record every pinned tool the expanded detection uses, with its version.

Writes refs/expanded/tools.lock.json: for each tool, where the binary is,
how it was installed, the version it reports and the version this project
pins. A mismatch is printed and the exit code is non-zero, so `make
tools-expanded` fails loudly instead of running a sample on a drifted tool.
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: name -> (path or command, version argv, pinned version substring, how installed)
TOOLS: dict[str, tuple[str, list[str], str, str]] = {
    "metaphlan (scoring lane)": (".venv-mpa3/bin/metaphlan", ["--version"], "3.1.0", "pip metaphlan==3.1.0 in .venv-mpa3"),
    "metaphlan (Jun23 baseline)": (".venv-mpa4/bin/metaphlan", ["--version"], "4.1.1", "pip metaphlan==4.1.1 in .venv-mpa4"),
    "metaphlan (Jan26 lane)": (".venv-mpa42/bin/metaphlan", ["--version"], "4.2.5",
                                "GitHub tag 4.2.6 (commit c55b299) in .venv-mpa42; the tag's setup.py still says 4.2.5"),
    "sample2markers": (".venv-mpa42/bin/sample2markers.py", ["--version"], "", "ships with MetaPhlAn 4.2.6"),
    "strainphlan": (".venv-mpa42/bin/strainphlan", ["--version"], "", "ships with MetaPhlAn 4.2.6"),
    "sylph": ("vendor/sylph/bin/sylph", ["--version"], "1.0.0", "built from source, vendor/sylph"),
    "motus": (".venv-motus/bin/motus", ["--version"], "4.1.0", "GitHub tag 4.1.0 in .venv-motus (python 3.12)"),
    "kraken2": ("kraken2", ["--version"], "2.1", "brew install kraken2"),
    "bracken": ("bracken", ["-v"], "", "brew install bracken"),
    "skani": ("skani", ["--version"], "0.3", "brew install skani"),
    "seqkit": ("seqkit", ["version"], "2.", "brew install seqkit"),
    "aria2c": ("aria2c", ["--version"], "1.3", "brew install aria2"),
    "fastANI": ("fastANI", ["--version"], "1.3", "brew install fastani"),
    "bwa": ("bwa", [], "0.7.19", "brew install bwa"),
    "bowtie2": ("bowtie2", ["--version"], "2.5", "brew install bowtie2"),
    "samtools": ("samtools", ["--version"], "1.", "brew install samtools"),
    "diamond": ("diamond", ["--version"], "2.", "brew install diamond"),
    "singlem": (".venv-singlem/bin/singlem", ["--version"], "0.21.4", "pip singlem==0.21.4 in .venv-singlem (python 3.12)"),
    "orfm": ("orfm", ["--version"], "0.7", "brew install brewsci/bio/orfm"),
    "smafa": (str(Path.home() / ".cargo/bin/smafa"), ["--version"], "", "cargo install smafa"),
    "mfqe": (str(Path.home() / ".cargo/bin/mfqe"), ["--version"], "", "cargo install mfqe"),
    "hmmsearch": ("hmmsearch", ["-h"], "HMMER 3", "brew install hmmer"),
    "spades": ("vendor/SPAdes-4.3.0-Darwin/bin/spades.py", ["--version"], "4.3.0", "GitHub release v4.3.0 Darwin x86_64 tarball, vendor/"),
}


def _version(path: str, argv: list[str]) -> str:
    exe = str(REPO / path) if not Path(path).is_absolute() and "/" in path else (shutil.which(path) or path)
    if not Path(exe).exists():
        return ""
    try:
        out = subprocess.run([exe, *argv], capture_output=True, text=True, check=False, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    import re
    text = (out.stdout or "") + (out.stderr or "")
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    # First a line that names a version, then any line carrying an x.y number.
    for ln in lines:
        if "version" in ln.lower() and re.search(r"\d+\.\d+", ln):
            return ln[:80]
    for ln in lines:
        if re.search(r"(^|\s|v)\d+\.\d+", ln) and not ln.startswith(("-", "CFLAGS", "Usage")):
            return ln[:80]
    return lines[0][:80] if lines else ""


def main() -> int:
    record: dict[str, dict[str, str | bool]] = {}
    bad = 0
    for name, (path, argv, pinned, how) in TOOLS.items():
        exe = str(REPO / path) if ("/" in path and not Path(path).is_absolute()) else (shutil.which(path) or path)
        present = Path(exe).exists()
        version = _version(path, argv) if present else ""
        ok = present and (pinned in version if pinned else True)
        if not ok:
            bad += 1
        record[name] = {"path": exe, "present": present, "reported_version": version, "pinned": pinned, "installed_by": how, "ok": ok}
        print(f"  {'ok ' if ok else 'BAD'} {name:28} {version[:50]!s:50} {'' if ok else '(missing or off pin)'}")
    out = REPO / "refs" / "expanded" / "tools.lock.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"locked_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
                               "platform": sys.platform, "tools": record}, indent=2) + "\n")
    print(f"\nwrote {out}; {len(TOOLS) - bad} of {len(TOOLS)} tools on pin")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
