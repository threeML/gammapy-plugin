#!/usr/bin/env bash
set -eo pipefail

OWNER="threeML"
REPO="gammapy-plugin"
SHA="${READTHEDOCS_GIT_COMMIT_HASH}"
TOKEN="${GITHUB_TOKEN}"

if [ -z "$TOKEN" ]; then
	echo "This is a Pull Request Build - cannot download artifacts safely"
	exit 0
fi

API="https://api.github.com/repos/$OWNER/$REPO"

echo "Searching artifacts for SHA: $SHA"

ARTIFACT_ID_ARRAY=$(
python3 - <<PY
import os, json, urllib.request

owner = "$OWNER"
repo = "$REPO"
sha = "$SHA"
token = "$TOKEN"

url = f"https://api.github.com/repos/{owner}/{repo}/actions/artifacts?per_page=100"

req = urllib.request.Request(url)
req.add_header("Accept", "application/vnd.github+json")
req.add_header("Authorization", f"Bearer {token}")

with urllib.request.urlopen(req) as r:
    data = json.load(r)

for a in data.get("artifacts", []):
    if a.get("expired"):
        continue
    name = a.get("name", "")
    if sha in name:
        print(a["id"])
	print(";")
PY
)

if [[ -z "$ARTIFACT_ID_ARRAY" ]]; then
  echo "No artifact found for SHA: $SHA"
  exit 1
fi
IFS=';' read -ra ARTIFACT_IDS <<< "$ARTIFACT_ID_ARRAY"
for ARTIFACT_ID in "${ARTIFACT_IDS[@]}"; do

	echo "Artifact ID: $ARTIFACT_ID"

	echo "Downloading artifact..."

	curl -fL \
	  -H "Authorization: Bearer $TOKEN" \
	  -H "Accept: application/vnd.github+json" \
	  "$API/actions/artifacts/$ARTIFACT_ID/zip" \
	  -o artifact.zip

	echo "Downloaded artifact.zip"
	if [[ $ARTIFACT_ID == *"notebooks"* ]];
	then
		echo "Found a notebooks artifact"
		mkdir -p docs/notebooks
		unzip -o artifact.zip -d docs/notebooks
	else
		echo "Found API stubs"

		mkdir -p docs/api
		unzip -o artifact.zip -d docs/api

	fi;

	rm artifact.zip
done

echo "Done unzipping"
