"""Add a configured Cloudflare Web Analytics beacon to built HTML pages."""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path


def insert(directory: Path, token: str) -> int:
    if not token:
        print("Web analytics is not configured; publishing pages without a beacon.")
        return 0
    if not re.fullmatch(r"[a-fA-F0-9]{32}", token):
        raise ValueError("Cloudflare Web Analytics token must be 32 hexadecimal characters")
    pages = sorted(directory.glob("*.html"))
    if not pages:
        raise ValueError(f"{directory}: no HTML pages found")
    snippet = (
        '<script type="module" src="https://static.cloudflareinsights.com/beacon.min.js" '
        f'data-cf-beacon=\'{{"token":"{token}"}}\'></script>'
    )
    for page in pages:
        html = page.read_text(encoding="utf-8")
        if html.count("</body>") != 1 or "static.cloudflareinsights.com/beacon.min.js" in html:
            raise ValueError(f"{page}: missing body close or analytics beacon already present")
        page.write_text(html.replace("</body>", f"  {snippet}\n</body>"), encoding="utf-8", newline="\n")
    print(f"Added Web Analytics to {len(pages)} HTML pages.")
    return len(pages)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/insert_site_analytics.py BUILD_DIRECTORY")
    insert(Path(sys.argv[1]), os.environ.get("CLOUDFLARE_WEB_ANALYTICS_TOKEN", ""))
