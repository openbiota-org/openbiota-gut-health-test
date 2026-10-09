"""MkDocs hook: every external link in page content opens in a new tab.

Markdown gives no way to say `target="_blank"` on a link, and the docs link
out to papers, licences and tools; none of those should replace the docs
tab. Runs on the rendered HTML of each page, so the built site is static
and needs no script for this. Relative links (the landing page, other docs
pages, anchors) are left alone.
"""

from __future__ import annotations

import re

# An <a ...> whose href is absolute http(s) and which has no target yet.
_EXTERNAL = re.compile(r'<a\s+(?![^>]*\btarget=)(?P<attrs>[^>]*\bhref="https?://[^"]+"[^>]*)>')


def on_page_content(html: str, page, config, files) -> str:  # noqa: ARG001 - MkDocs hook signature
    return _EXTERNAL.sub(r'<a target="_blank" rel="noopener noreferrer" \g<attrs>>', html)
