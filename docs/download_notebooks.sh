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
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


owner = os.environ["OWNER"]
repo = os.environ["REPO"]
sha = os.environ["SHA"]
token = os.environ["TOKEN"]

api = f"https://api.github.com/repos/{owner}/{repo}"

github_headers = {
    "Accept": "application/vnd.github+json",
    "Authorization": f"Bearer {token}",
    "X-GitHub-Api-Version": "2022-11-28",
}


def github_request(url):
    """Make an authenticated request to the GitHub API."""
    request = urllib.request.Request(
        url,
        headers=github_headers,
    )
    return urllib.request.urlopen(request)


class NoRedirect(urllib.request.HTTPRedirectHandler):
    """Prevent urllib from automatically following redirects."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


no_redirect_opener = urllib.request.build_opener(NoRedirect)


def get_artifact_download_url(url):
    """
    Ask GitHub for the artifact download URL without following
    the redirect.

    The GitHub Authorization header is sent only to api.github.com.
    """
    request = urllib.request.Request(
        url,
        headers=github_headers,
    )

    try:
        no_redirect_opener.open(request)

    except urllib.error.HTTPError as exc:
        if exc.code != 302:
            raise

        location = exc.headers.get("Location")

        if not location:
            raise RuntimeError(
                "GitHub returned HTTP 302 without a Location header"
            )

        return location

    raise RuntimeError(
        "Expected GitHub artifact endpoint to return HTTP 302"
    )


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
    print("Requesting artifact download URL...")

    github_download_url = (
        f"{api}/actions/artifacts/{artifact_id}/zip"
    )

    # GitHub returns a short-lived signed URL via HTTP 302.
    signed_download_url = get_artifact_download_url(
        github_download_url
    )

    print("Downloading artifact...")

    zip_path = Path(f"artifact-{artifact_id}.zip")

    try:
        # IMPORTANT:
        # Do NOT send the GitHub Authorization header here.
        # The signed URL already contains its own authentication.
        with urllib.request.urlopen(signed_download_url) as response:
            with zip_path.open("wb") as output:
                shutil.copyfileobj(response, output)

        print(f"Downloaded {zip_path}")

        if "notebooks" in artifact_name.lower():
            print("Found a notebooks artifact")
            output_dir = Path("docs/notebooks")
        else:
            print("Found API stubs")
            output_dir = Path("docs/api")

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        print(f"Extracting to {output_dir}")

        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(output_dir)

        print(f"Extracted {artifact_name}")

    finally:
        zip_path.unlink(missing_ok=True)


print("Artifact download complete")
PY
