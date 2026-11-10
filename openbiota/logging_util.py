"""Progress reporting to stderr plus a persistent run log.

The build spec calls a 30-minute silent run "a bug report waiting to happen",
so every stage emits a timestamped line and all of them are accumulated for
``run.log``. What reaches the terminal comes in two sizes:

* **detail** (``Reporter(detail=True)``, `openbiota run --verbose`, and every
  maintenance command): every line, timestamped, as the log has it.
* **short** (``detail=False``, the default for `openbiota run`): one line per
  stage that took a moment, with how long it took; warnings and failures in
  one line each. On a terminal the stage in progress is a single line that
  redraws once a second with the elapsed time and the latest detail, so a
  twenty-minute stage is visibly alive without scrolling anything; work
  between stages shows the same way, as "working". Where a line cannot be
  redrawn (a pipe, a log file), a stage that is still going after a couple
  of seconds is announced once, so the last line of a log always says what
  is running.
"""

from __future__ import annotations

import os
import re
import shutil
import sys
import threading
import time
import weakref
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Final

_ANSI = re.compile(r"\033\[[0-9;]*m")

#: In the short form, a stage quicker than this leaves no line.
STAGE_SHOWN_AFTER_S = 1.0
#: Where the running stage cannot be a redrawn line, it is announced after this long.
STAGE_ANNOUNCED_AFTER_S = 2.0
#: On a terminal, work between two stages gets a status line after this long.
IDLE_SHOWN_AFTER_S = 2.0
#: The column the duration is aligned to in the short form's stage lines.
_STAGE_NAME_WIDTH = 58


def _supports_colour(stream: object) -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    isatty = getattr(stream, "isatty", None)
    return bool(isatty and isatty())


def _redraws(stream: object) -> bool:
    """Whether a line on this stream can be redrawn in place."""
    isatty = getattr(stream, "isatty", None)
    return bool(isatty and isatty()) and os.environ.get("TERM", "") != "dumb"


def _first_line(text: str) -> str:
    for line in _ANSI.sub("", text or "").splitlines():
        if line.strip():
            return line.strip()
    return ""


_NOT_FOUND = re.compile(r"(?:^|[\s:'\"`(])([A-Za-z0-9_.+-]{2,40}): (?:command )?not found")
_NO_FILE = re.compile(r"No such file or directory: '([^']+)'")
_MISSING_MODULE = re.compile(r"No module named '([^']+)'")
#: How a process ended by the kernel shows up: the shell's "Killed" (macOS adds
#: ": 9") as the last word, a SIGKILL in a Python error, exit status 137.
_KILLED = re.compile(r"(?:^|[\s:(])Killed(?:: 9)?\s*$|\b\d+ Killed\b|SIGKILL|died with signal 9\b|"
                     r"exit(?: status)? 137\b|exit status -9\b", re.MULTILINE)
OUT_OF_MEMORY: Final = ("stopped by the operating system - almost always because the machine ran out of "
                        "memory (docs/INSTALL.md lists what each lane needs)")
_EXCEPTION_PREFIX = re.compile(r"^(?:(?:failed|error)\s*:\s*)?(?:[A-Za-z_.]*(?:Error|Exception)\s*:\s*)+", re.IGNORECASE)


def plain_reason(detail: str, *, limit: int = 160) -> str:
    """One readable line from whatever a failed tool left behind.

    Tools fail with a traceback wrapped around a subprocess error wrapped
    around a shell message. The part a person can act on is usually the
    innermost one - which program was not found - so that is looked for
    first; otherwise the first line that says something, without the
    exception class names in front of it.
    """
    text = _ANSI.sub("", detail or "").replace("\\n", "\n")
    found = _NOT_FOUND.search(text)
    if found:
        return f"the program '{found.group(1)}' is not installed - run ./install.sh to add it"
    missing = _NO_FILE.search(text)
    if missing:
        return f"missing file or program: {missing.group(1)} - run ./install.sh"
    module = _MISSING_MODULE.search(text)
    if module:
        return f"a Python package is missing ({module.group(1)}) - run ./install.sh"
    if _KILLED.search(text):
        return OUT_OF_MEMORY
    for line in text.splitlines():
        line = line.strip().strip("b'\"")
        if not line or line.startswith(("Traceback", "File ", "STDERR was", "During handling")):
            continue
        line = _EXCEPTION_PREFIX.sub("", line).strip() or line
        return line[:limit] + ("..." if len(line) > limit else "")
    return "failed without a message; see run.log"


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


def _tick(ref: weakref.ReferenceType[Reporter]) -> None:
    """The once-a-second refresh. Holds no reference: it ends with its reporter."""
    while True:
        time.sleep(1.0)
        reporter = ref()
        if reporter is None:
            return
        with reporter._lock:  # noqa: SLF001
            if reporter._closed or (reporter._live is None and not reporter._redraw):  # noqa: SLF001
                reporter._ticker = None  # noqa: SLF001
                return
            reporter._refresh()  # noqa: SLF001
        del reporter


@dataclass(slots=True, weakref_slot=True)
class Reporter:
    """Stage-aware progress reporter.

    Lines go to stderr as they happen (so a long run is never silent) and
    are buffered so the CLI can persist them to ``results/<sample>/run.log``.
    ``verbose=False`` prints nothing; ``detail=False`` prints the short form
    described at the top of this module. The log gets every line either way.
    """

    verbose: bool = True
    detail: bool = True
    _t0: float = field(default_factory=time.monotonic)
    _lines: list[str] = field(default_factory=list)
    _style: Style = field(default_factory=lambda: Style(_supports_colour(sys.stderr)))
    _lock: threading.RLock = field(default_factory=threading.RLock)
    #: Short form only: whether the running stage is one redrawn line.
    _redraw: bool = field(default_factory=lambda: _redraws(sys.stderr))
    _depth: int = 0
    _live: tuple[str, float] | None = None
    _live_detail: str = ""
    _live_drawn: bool = False
    _announced: bool = False
    _idle_since: float = 0.0
    _closed: bool = False
    _ticker: threading.Thread | None = None
    _last_failure: str = ""
    _warned: set[str] = field(default_factory=set)

    @property
    def style(self) -> Style:
        return self._style

    @property
    def elapsed(self) -> float:
        return time.monotonic() - self._t0

    # -- the log ----------------------------------------------------------- #

    def _log(self, text: str) -> str:
        """Append to the log (plain text, timestamped); returns the stamp."""
        stamp = f"[{self.elapsed:7.1f}s]"
        self._lines.append(f"{stamp} {_ANSI.sub('', text)}")
        return stamp

    def _emit(self, text: str, *, to_stderr: bool = True) -> None:
        with self._lock:
            stamp = self._log(text)
            if not (to_stderr and self.verbose):
                return
            if self.detail:
                print(f"{self._style.dim(stamp)} {text}", file=sys.stderr, flush=True)
            else:
                # Progress detail: beside the running stage, never a line of its own.
                self._live_detail = _first_line(text)

    # -- the terminal, short form --------------------------------------------- #

    def _line(self, text: str) -> None:
        """A permanent line; the redrawn stage line steps aside and comes back."""
        self._clear_live()
        print(text, file=sys.stderr, flush=True)
        self._draw_live()

    def _clear_live(self) -> None:
        if self._live_drawn:
            sys.stderr.write("\r\033[2K")
            sys.stderr.flush()
            self._live_drawn = False

    def _draw_live(self, idle: tuple[str, float] | None = None) -> None:
        current = self._live or idle
        if current is None or not self._redraw or not self.verbose or self._closed:
            return
        name, t0 = current
        seconds = int(time.monotonic() - t0)
        took = f"{seconds}s" if seconds < 60 else _human_duration(seconds)
        room = max(20, shutil.get_terminal_size((100, 24)).columns - 1)
        name = name[: max(8, room - len(took) - 6)]
        used = 4 + len(name) + 2 + len(took)
        latest = self._live_detail[: room - used - 2] if room - used - 2 >= 12 else ""
        sys.stderr.write(
            "\r\033[2K  " + self._style.cyan("·") + f" {name}  " + self._style.dim(took)
            + ("  " + self._style.dim(latest) if latest else "")
        )
        sys.stderr.flush()
        self._live_drawn = True

    def _refresh(self) -> None:
        """Once a second: redraw the running stage's line, or announce it once."""
        if not self.verbose or self._closed:
            return
        now = time.monotonic()
        if self._live is None:
            # Between stages. On a terminal that is worth a line too, or a
            # minute of unstaged work looks like a hang.
            if self._redraw and self._idle_since and now - self._idle_since >= IDLE_SHOWN_AFTER_S:
                self._draw_live(("working", self._idle_since))
        elif self._redraw:
            self._draw_live()
        elif not self._announced and now - self._live[1] >= STAGE_ANNOUNCED_AFTER_S:
            print(f"  · {self._live[0]} …", file=sys.stderr, flush=True)
            self._announced = True

    def close(self) -> None:
        """The run's work is over: wipe the status line and stop refreshing it."""
        with self._lock:
            self._clear_live()
            self._live = None
            self._closed = True

    def _begin(self, name: str) -> bool:
        """Open a stage in the short form; True when it is a top-level one."""
        with self._lock:
            self._log(f"· {name} …")
            top = self._depth == 0
            self._depth += 1
            if not top:
                self._live_detail = name   # a step inside the running stage: beside it, not under it
            if top and self.verbose and not self._closed:
                self._live, self._live_detail, self._announced = (name, time.monotonic()), "", False
                self._draw_live()
                if self._ticker is None:
                    self._ticker = threading.Thread(target=_tick, args=(weakref.ref(self),),
                                                    name="openbiota-status", daemon=True)
                    self._ticker.start()
            return top

    def _finish(self, name: str, t0: float, top: bool, failed: BaseException | None) -> None:
        took = time.monotonic() - t0
        with self._lock:
            self._depth = max(0, self._depth - 1)
            if failed is None:
                self._log(f"✓ {name} done in {_human_duration(took)}")
            else:
                self._log(f"✗ {name} failed after {took:.1f}s: {type(failed).__name__}: {failed}")
            if not top:
                return
            self._clear_live()
            self._live, self._live_detail = None, ""
            self._idle_since = time.monotonic()
            if not self.verbose:
                return
            if failed is None:
                if took >= STAGE_SHOWN_AFTER_S or self._announced:
                    print(f"  {self._style.green('✓')} {name:<{_STAGE_NAME_WIDTH}} {_human_duration(took):>8}",
                          file=sys.stderr, flush=True)
            elif isinstance(failed, Exception):
                # Not for Ctrl-C or SystemExit: whoever handles those says so.
                message = _first_line(str(failed))
                self._last_failure = message[:60]
                print(f"  {self._style.red('✗')} {name}: {plain_reason(f'{type(failed).__name__}: {failed}')}",
                      file=sys.stderr, flush=True)

    # -- what callers say --------------------------------------------------- #

    def info(self, text: str) -> None:
        self._emit(text)

    def step(self, text: str) -> None:
        self._emit(self._style.cyan("· ") + text)

    def ok(self, text: str) -> None:
        self._emit(self._style.green("✓ ") + text)

    def warn(self, text: str) -> None:
        if self.detail:
            self._emit(self._style.yellow("! ") + text)
            return
        with self._lock:
            self._log("! " + text)
            first = _first_line(text)
            # A handler explaining the failure the stage line just reported
            # adds nothing; nor does the same warning a second time.
            explained = bool(self._last_failure) and self._last_failure in _ANSI.sub("", text)
            self._last_failure = ""
            if not self.verbose or explained or first in self._warned:
                return
            self._warned.add(first)
            self._line(f"  {self._style.yellow('!')} {first[:220]}{'...' if len(first) > 220 else ''}")

    def error(self, text: str) -> None:
        if self.detail:
            self._emit(self._style.red("✗ ") + text)
            return
        with self._lock:
            self._log("✗ " + text)
            if self.verbose:
                self._line(f"  {self._style.red('✗')} {_first_line(text)}")

    def headline(self, text: str) -> None:
        """A line the short form shows as written; the detailed form stamps it like any other."""
        if self.detail:
            self._emit(text)
            return
        with self._lock:
            self._log(text)
            if self.verbose:
                self._line(text)

    def record(self, text: str) -> None:
        """Add to run.log without printing to stderr."""
        self._emit(text, to_stderr=False)

    @contextmanager
    def stage(self, name: str) -> Iterator[None]:
        if self.detail:
            self.step(f"{name} …")
            t0 = time.monotonic()
            try:
                yield
            except BaseException:
                self.error(f"{name} failed after {time.monotonic() - t0:.1f}s")
                raise
            self.ok(f"{name} done in {_human_duration(time.monotonic() - t0)}")
            return
        top = self._begin(name)
        t0 = time.monotonic()
        try:
            yield
        except BaseException as exc:
            self._finish(name, t0, top, exc)
            raise
        self._finish(name, t0, top, None)

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
