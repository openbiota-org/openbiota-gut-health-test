"""MkDocs hook: every path written into the site is lower-case.

MkDocs names each output file after its Markdown source, and the sources
follow the repository convention of upper-case names (`QUICKSTART.md`), which
would give URLs like `/docs/QUICKSTART.html`. URLs should not shout, so the
destination of every file - pages, folders, images, theme assets - is
lower-cased here before anything is rendered. MkDocs builds the navigation,
internal links and sitemap from `File.url`, so they all follow automatically;
the sources keep their names and the "edit this page" links still point at
them.

A post-build check walks the finished site and, because the docs build runs
with `strict: true`, fails the build if anything upper-case slipped through.
"""

from __future__ import annotations

import logging
import os

log = logging.getLogger("mkdocs.hooks.lowercase_urls")


def on_files(files, config):  # noqa: ARG001 - MkDocs hook signature
    for file in files:
        lowered = file.dest_uri.lower()
        if lowered != file.dest_uri:
            file.dest_uri = lowered
            # `url` and `abs_dest_path` are cached from `dest_uri`; drop them
            # so they are recomputed from the lower-cased path.
            file.__dict__.pop("url", None)
            file.__dict__.pop("abs_dest_path", None)
    return files


def on_post_build(config) -> None:
    site_dir = config["site_dir"]
    offenders = []
    for root, dirs, names in os.walk(site_dir):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for name in dirs + [n for n in names if not n.startswith(".")]:
            rel = os.path.relpath(os.path.join(root, name), site_dir)
            if rel != rel.lower():
                offenders.append(rel)
    for rel in sorted(offenders):
        log.warning("Upper-case path in the built site: %s", rel)
