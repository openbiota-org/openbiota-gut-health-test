"""Minimal, dependency-free HTTP with retry and backoff.

Only used for the UniProt REST fetches, which happen once and are then cached
on disk. Kept in its own module so the retry policy is testable in isolation.
"""

from __future__ import annotations

import http.client
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from openbiota import __version__
from openbiota.errors import ReferenceFetchError

USER_AGENT: Final = f"openbiota/{__version__} (OpenBiota Gut Health Test; +https://openbiota.com)"

#: HTTP statuses worth retrying: transient server faults and rate limiting.
RETRY_STATUSES: Final = frozenset({408, 425, 429, 500, 502, 503, 504})


@dataclass(frozen=True, slots=True)
class Response:
    url: str
    status: int
    body: bytes
    headers: dict[str, str]

    def header(self, name: str) -> str | None:
        lowered = name.lower()
        for key, value in self.headers.items():
            if key.lower() == lowered:
                return value
        return None

    def text(self) -> str:
        return self.body.decode("utf-8", errors="replace")


def get(
    url: str,
    *,
    timeout: float = 120.0,
    attempts: int = 5,
    backoff: float = 2.0,
) -> Response:
    """GET ``url``, retrying transient failures with exponential backoff."""
    last_error: str = "no attempt made"
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return Response(
                    url=url,
                    status=int(response.status),
                    body=response.read(),
                    headers=dict(response.headers.items()),
                )
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code} {exc.reason}"
            if exc.code not in RETRY_STATUSES:
                raise ReferenceFetchError(f"{url}: {last_error}") from exc
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"

        if attempt < attempts:
            delay = backoff ** (attempt - 1)
            time.sleep(delay)

    raise ReferenceFetchError(
        f"{url}: giving up after {attempts} attempts — {last_error}. "
        "Check network connectivity, or reuse a cached reference set by omitting --refresh-refs."
    )


def download_to(
    url: str, dest: Path, *, timeout: float = 600.0, attempts: int = 6, chunk: int = 1 << 22
) -> Path:
    """Stream ``url`` to ``dest``, resuming a partial download on retry.

    For multi-hundred-megabyte reference genomes, reading the whole body into
    memory and retrying from zero on a dropped connection is not acceptable.
    This writes to ``dest.part`` and sends a ``Range`` header on retry.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last_error = "no attempt made"
    for attempt in range(1, attempts + 1):
        have = part.stat().st_size if part.is_file() else 0
        headers = {"User-Agent": USER_AGENT}
        if have:
            headers["Range"] = f"bytes={have}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                if have and response.status != 206:
                    have = 0  # server ignored the range; start over
                mode = "ab" if have else "wb"
                with part.open(mode) as handle:
                    while True:
                        block = response.read(chunk)
                        if not block:
                            break
                        handle.write(block)
            part.replace(dest)
            return dest
        except urllib.error.HTTPError as exc:
            last_error = f"HTTP {exc.code} {exc.reason}"
            if exc.code == 416 and part.is_file():
                part.replace(dest)  # range not satisfiable: already complete
                return dest
            if exc.code not in RETRY_STATUSES:
                raise ReferenceFetchError(f"{url}: {last_error}") from exc
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
        if attempt < attempts:
            time.sleep(min(60.0, 2.0 ** (attempt - 1)))
    raise ReferenceFetchError(f"{url}: giving up after {attempts} attempts — {last_error}")


def build_url(base: str, params: dict[str, str | int]) -> str:
    return f"{base}?{urllib.parse.urlencode(params)}"


#: RFC 8288 Link header entries. Splitting the header on "," is wrong: the
#: UniProt URL itself contains commas in its `fields=` parameter.
_LINK_RE: Final = re.compile(r"<(?P<url>[^>]+)>\s*;\s*(?P<params>[^<]*)")


def next_page_url(response: Response) -> str | None:
    """Extract the cursor-paging ``rel="next"`` link from a UniProt response."""
    link = response.header("Link")
    if not link:
        return None
    for match in _LINK_RE.finditer(link):
        params = match.group("params").replace("'", '"').replace(" ", "")
        if 'rel="next"' in params or "rel=next" in params:
            return match.group("url").strip()
    return None
