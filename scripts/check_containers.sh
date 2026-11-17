#!/usr/bin/env bash
# Pull every biocontainers image listed in nextflow.config and report which tags fail.
# For a failed image, list tags that exist on quay.io for the same tool and version
# so you can paste the right one into nextflow.config.
#
# Usage (needs Docker running):  bash scripts/check_containers.sh
set -u
CONFIG="${1:-nextflow.config}"

grep -oE "quay.io/biocontainers/[^:' ]+:[^' ]+" "$CONFIG" | sort -u | while read -r img; do
  rest="${img#quay.io/biocontainers/}"
  tool="${rest%%:*}"
  tag="${rest#*:}"
  ver="${tag%%--*}"
  if docker pull "$img" >/dev/null 2>&1; then
    echo "OK      $img"
  else
    echo "FAILED  $img"
    echo "        tags on quay.io for '$tool' containing '$ver':"
    python3 - "$tool" "$ver" <<'PY' || echo "        (could not query quay.io)"
import json, sys, urllib.request
tool, ver = sys.argv[1], sys.argv[2]
found, page = [], 1
while page <= 20:
    url = f"https://quay.io/api/v1/repository/biocontainers/{tool}/tag/?limit=100&page={page}&onlyActiveTags=true"
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    found += [t["name"] for t in data.get("tags", []) if ver in t["name"]]
    if not data.get("has_additional"):
        break
    page += 1
for name in sorted(set(found)) or ["(none found - try a different version)"]:
    print("          quay.io/biocontainers/%s:%s" % (tool, name) if not name.startswith("(") else "          " + name)
PY
  fi
done
