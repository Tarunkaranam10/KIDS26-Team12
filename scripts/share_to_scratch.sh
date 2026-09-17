#!/bin/bash
# Stage the KIDS26 Team 12 public TCGA dataset to shared scratch for collaborators.
# Submitted to LSF via scripts/submit_share.sh; safe to rerun (rsync is incremental).
set -euo pipefail

SRC="/research/groups/shelagrp/home/esavage/KIDS26/KIDS26-Team12"
DEST="/lustre_scratch/shared_scratch/kids26_team12_share"

echo "host=$(hostname) start=$(date -Is)"
echo "src=$SRC"
echo "dest=$DEST"

mkdir -p "$DEST"/data/raw "$DEST"/data/processed "$DEST"/config "$DEST"/scripts "$DEST"/docs

# Bulk payload. --partial keeps resumable chunks; -h prints human-readable sizes.
rsync -av --partial --human-readable \
  --exclude '*.part' \
  "$SRC/data/raw/" "$DEST/data/raw/"

rsync -av --partial --human-readable \
  --exclude '*.part' \
  "$SRC/data/processed/" "$DEST/data/processed/"

# Metadata, manifests, code and documentation needed to interpret and regenerate.
rsync -av "$SRC/config/" "$DEST/config/"
rsync -av "$SRC/scripts/" "$DEST/scripts/" --exclude '__pycache__'
rsync -av "$SRC/docs/" "$DEST/docs/"
rsync -av "$SRC/README.md" "$DEST/repo_README.md"
rsync -av "$SRC/README_BIOHACKATHON.md" "$DEST/repo_README_BIOHACKATHON.md"

# Collaborator-facing dataset README.
rsync -av "$SRC/share/README.md" "$DEST/README.md"

echo "rsync complete=$(date -Is)"

# Independent verification against the sidecars written at download time.
# NOTE: some upstream source manifests are CRLF, and the trailing \r can leak into the
# filename field of a generated .sha256 sidecar. Strip CR before checking.
cd "$DEST/data/raw/pancanatlas"
echo "=== sha256 verification (data/raw/pancanatlas) ==="
cat ./*.sha256 | tr -d '\r' > /tmp/kids26_share_expected.$$.sha256
sha256sum -c /tmp/kids26_share_expected.$$.sha256
rm -f /tmp/kids26_share_expected.$$.sha256

echo "=== sha256 verification (data/raw/gdc, beta + masked-idat) ==="
cd "$DEST/data/raw/gdc"
find . -name '*.sha256' -print0 | while IFS= read -r -d '' s; do
  ( cd "$(dirname "$s")" && tr -d '\r' < "$(basename "$s")" | sha256sum -c - )
done

echo "=== processed beta.tsv sha256 ==="
sha256sum "$DEST/data/processed/beta.tsv"
echo "expected from beta.provenance.json:"
grep '"sha256"' "$DEST/data/processed/beta.provenance.json"

chmod -R g+rX,o+rX "$DEST" || true

echo "=== final size ==="
du -sh "$DEST"
echo "done=$(date -Is)"
