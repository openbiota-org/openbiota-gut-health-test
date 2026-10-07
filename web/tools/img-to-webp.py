#!/usr/bin/env python3
"""
Convert PNG/JPEG images to WebP, writing <name>.webp next to each source.

Usage:
  python3 tools/img-to-webp.py <file-or-dir> [more...] [--quality 82] [--max-width N] [--lossless] [--force]

Examples:
  python3 tools/img-to-webp.py wwwroot/img/portfolio            # all PNG/JPG in the folder
  python3 tools/img-to-webp.py wwwroot/img/portfolio/newco.png  # one new thumbnail
  python3 tools/img-to-webp.py wwwroot/img/hero.jpg --quality 78 --max-width 1920

Sources are never modified or deleted. Existing .webp files are skipped unless the
source is newer or --force is given. Requires Pillow with WebP support (python3 -m pip install pillow).
"""
import argparse
import sys
from pathlib import Path

try:
    from PIL import Image, features
except ImportError:
    sys.exit("Pillow is required: python3 -m pip install pillow")

if not features.check("webp"):
    sys.exit("This Pillow build lacks WebP support; reinstall Pillow (python3 -m pip install --upgrade pillow).")

SOURCE_EXTS = {".png", ".jpg", ".jpeg"}


def collect(paths):
    for raw in paths:
        p = Path(raw)
        if p.is_dir():
            yield from sorted(f for f in p.iterdir() if f.suffix.lower() in SOURCE_EXTS)
        elif p.is_file() and p.suffix.lower() in SOURCE_EXTS:
            yield p
        else:
            print(f"skip (not a PNG/JPEG or not found): {raw}", file=sys.stderr)


def convert(src, quality, max_width, lossless, force):
    """Write src as .webp; return (bytes_before, bytes_after) or None if skipped."""
    dst = src.with_suffix(".webp")
    if dst.exists() and not force and dst.stat().st_mtime >= src.stat().st_mtime:
        print(f"up to date  {dst}")
        return None

    with Image.open(src) as im:
        # Preserve transparency; flatten palette images so WebP gets real pixels.
        if im.mode in ("P", "LA"):
            im = im.convert("RGBA")
        elif im.mode not in ("RGB", "RGBA"):
            im = im.convert("RGB")

        if max_width and im.width > max_width:
            ratio = max_width / im.width
            im = im.resize((max_width, round(im.height * ratio)), Image.LANCZOS)

        save_kwargs = {"method": 6}  # slowest / best compression
        if lossless:
            save_kwargs["lossless"] = True
        else:
            save_kwargs["quality"] = quality
        im.save(dst, "WEBP", **save_kwargs)

    before, after = src.stat().st_size, dst.stat().st_size
    pct = (1 - after / before) * 100 if before else 0
    print(f"{src.name:40s} {before/1024:7.1f} KB -> {after/1024:7.1f} KB  ({pct:5.1f}% smaller)  {dst}")
    return before, after


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="image files and/or directories")
    ap.add_argument("--quality", type=int, default=82, help="lossy quality 1-100 (default 82)")
    ap.add_argument("--max-width", type=int, default=0, help="downscale wider images to this width (default: keep size)")
    ap.add_argument("--lossless", action="store_true", help="lossless WebP (good for flat-color logos)")
    ap.add_argument("--force", action="store_true", help="re-encode even if the .webp is already up to date")
    args = ap.parse_args()

    total_before = total_after = converted = 0
    for src in collect(args.paths):
        result = convert(src, args.quality, args.max_width, args.lossless, args.force)
        if result:
            converted += 1
            total_before += result[0]
            total_after += result[1]

    if converted:
        pct = (1 - total_after / total_before) * 100 if total_before else 0
        print(f"\n{converted} file(s): {total_before/1024:.0f} KB -> {total_after/1024:.0f} KB ({pct:.0f}% smaller)")
    else:
        print("\nNothing to convert.")


if __name__ == "__main__":
    main()
