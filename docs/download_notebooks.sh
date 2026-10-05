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

export OWNER
export REPO
export SHA
export TOKEN

python3 - <<'PY'
import json
import os
import shutil
import urllib.request
import zipfile
from pathlib import Path


owner = os.environ["OWNER"]
repo = os.environ["REPO"]
sha = os.environ["SHA"]
token = os.environ["TOKEN"]

api = f"https://api.github.com/repos/{owner}/{repo}"

headers = {
    "Accept": "application/vnd.github+json",
    "Authorization": f"Bearer {token}",
    "X-GitHub-Api-Version": "2022-11-28",
}


def github_request(url):
    request = urllib.request.Request(url, headers=headers)
    return urllib.request.urlopen(request)


print(f"Searching artifacts for SHA: {sha}")

url = f"{api}/actions/artifacts?per_page=100"

with github_request(url) as response:
    data = json.load(response)


artifacts = []

for artifact in data.get("artifacts", []):
    if artifact.get("expired"):
        continue

    name = artifact.get("name", "")

    if sha in name:
        artifacts.append(artifact)


if not artifacts:
    print(f"No artifact found for SHA: {sha}")
    raise SystemExit(1)


print(f"Found {len(artifacts)} artifact(s)")


for artifact in artifacts:
    artifact_id = artifact["id"]
    artifact_name = artifact.get("name", "")

    print(f"Artifact ID: {artifact_id}")
    print(f"Artifact name: {artifact_name}")
    print("Downloading artifact...")

    zip_path = Path(f"artifact-{artifact_id}.zip")
    download_url = f"{api}/actions/artifacts/{artifact_id}/zip"

    try:
        with github_request(download_url) as response:
            with zip_path.open("wb") as output:
                shutil.copyfileobj(response, output)

        print(f"Downloaded {zip_path}")

        if "notebooks" in artifact_name.lower():
            print("Found a notebooks artifact")
            output_dir = Path("docs/notebooks")
        else:
            print("Found API stubs")
            output_dir = Path("docs/api")

        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Extracting to {output_dir}")

        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(output_dir)

        print(f"Extracted {artifact_name}")

    finally:
        zip_path.unlink(missing_ok=True)


print("Artifact download complete")
PY
