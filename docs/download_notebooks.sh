#!/usr/bin/env bash
set -euo pipefail

OWNER="threeML"
REPO="gammapy-plugin"
SHA="${READTHEDOCS_GIT_COMMIT_HASH}"
TOKEN="${GITHUB_TOKEN}"

API="https://api.github.com/repos/$OWNER/$REPO"

# Find artifact ID whose name contains the SHA
#

ARTIFACT_ID=$(
  curl -fsSL \
    -H "Authorization: Bearer $TOKEN" \
    -H "Accept: application/vnd.github+json" \
    "$API/actions/artifacts?per_page=100" \
  | jq -r --arg SHA "$SHA" '
      .artifacts[]
      | select(.expired == false)
      | select(.name | contains($SHA))
      | .id
    ' \
  | head -n1
)

if [[ -z "$ARTIFACT_ID" ]]; then
  echo "No artifact found for SHA: $SHA"
  exit 1
fi

echo "Artifact ID: $ARTIFACT_ID"

# Download artifact ZIP
curl -fL \
  -H "Authorization: Bearer $TOKEN" \
  -H "Accept: application/vnd.github+json" \
  "$API/actions/artifacts/$ARTIFACT_ID/zip" \
  -o artifact.zip

echo "Downloaded artifact.zip"
mkdir -p docs/notebooks
unzip artifact.zip -d docs/notebooks
rm artifac.zip
echo "Done unzipping"
