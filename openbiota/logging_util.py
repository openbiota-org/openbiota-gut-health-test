"""Progress reporting to stderr plus a persistent run log.

The build spec calls a 30-minute silent run "a bug report waiting to happen",
so every stage emits a timestamped line to stderr and the same line is
accumulated for ``run.log``.
"""

from __future__ import annotations

import os
import sys
import threading
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path


def _supports_colour(stream: object) -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    isatty = getattr(stream, "isatty", None)
    return bool(isatty and isatty())


@dataclass(slots=True)
class Style:
    """ANSI styling that degrades to plain text when stderr is not a TTY."""

    enabled: bool = False

    def _wrap(self, code: str, text: str) -> str:
        return f"\033[{code}m{text}\033[0m" if self.enabled else text

    def bold(self, text: str) -> str:
        return self._wrap("1", text)

    def dim(self, text: str) -> str:
        return self._wrap("2", text)

    def red(self, text: str) -> str:
        return self._wrap("31", text)

    def green(self, text: str) -> str:
        return self._wrap("32", text)

    def yellow(self, text: str) -> str:
        return self._wrap("33", text)

    def cyan(self, text: str) -> str:
        return self._wrap("36", text)


@dataclass(slots=True)
class Reporter:
    """Stage-aware progress reporter.

    Lines go to stderr immediately (so a long run is never silent) and are
    buffered so the CLI can persist them to ``results/<sample>/run.log``.
    """

    verbose: bool = True
    _t0: float = field(default_factory=time.monotonic)
    _lines: list[str] = field(default_factory=list)
    _style: Style = field(default_factory=lambda: Style(_supports_colour(sys.stderr)))
    _lock: threading.Lock = field(default_factory=threading.Lock)

    @property
    def style(self) -> Style:
        return self._style

    @property
    def elapsed(self) -> float:
        return time.monotonic() - self._t0

    def _emit(self, text: str, *, to_stderr: bool = True) -> None:
        stamp = f"[{self.elapsed:7.1f}s]"
        line = f"{stamp} {text}"
        with self._lock:
            self._lines.append(line)
            if to_stderr and self.verbose:
                print(f"{self._style.dim(stamp)} {text}", file=sys.stderr, flush=True)

    def info(self, text: str) -> None:
        self._emit(text)

    def step(self, text: str) -> None:
        self._emit(self._style.cyan("· ") + text)

    def ok(self, text: str) -> None:
        self._emit(self._style.green("✓ ") + text)

    def warn(self, text: str) -> None:
        self._emit(self._style.yellow("! ") + text)

    def error(self, text: str) -> None:
        self._emit(self._style.red("✗ ") + text)

    def record(self, text: str) -> None:
        """Add to run.log without printing to stderr."""
        self._emit(text, to_stderr=False)

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        self.step(f"{name} …")
        t0 = time.monotonic()
        try:
            yield
        except BaseException:
            self.error(f"{name} failed after {time.monotonic() - t0:.1f}s")
            raise
        self.ok(f"{name} done in {_human_duration(time.monotonic() - t0)}")

    @contextmanager
    def heartbeat(self, label: str, probe: Callable[[], str], interval: float = 20.0) -> Iterator[None]:
        """Emit ``label`` plus ``probe()`` every ``interval`` seconds.

        Used to keep the long DIAMOND stages visibly alive. ``probe`` must be
        cheap and must not raise; it typically reports output-file growth.
        """
        stop = threading.Event()

        def loop() -> None:
            t0 = time.monotonic()
            while not stop.wait(interval):
                try:
                    detail = probe()
                except OSError:
                    detail = "n/a"
                self.info(f"  {label} running {_human_duration(time.monotonic() - t0)} — {detail}")

        thread = threading.Thread(target=loop, name="openbiota-heartbeat", daemon=True)
        thread.start()
        try:
            yield
        finally:
            stop.set()
            thread.join(timeout=1.0)

    def write_log(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with self._lock:
            path.write_text("\n".join(self._lines) + "\n", encoding="utf-8")


def _human_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, secs = divmod(int(seconds), 60)
    if minutes < 60:
        return f"{minutes}m{secs:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{secs:02d}s"


def human_duration(seconds: float) -> str:
    return _human_duration(seconds)


def human_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0 or unit == "TB":
            return f"{n:.1f}{unit}" if unit != "B" else f"{int(n)}B"
        n /= 1024.0
    return f"{n:.1f}TB"
