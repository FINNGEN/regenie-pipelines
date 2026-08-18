#!/bin/bash

usage() {
  echo "Usage: $0 --version <finngen_tag> [--image <name>] [--push] [--registry refinery|sandbox] [--base-regenie-docker <image>]"
  echo "  --version              required. Suffix appended to regenie's own VERSION file to build the pushed tag (e.g. cond_firth)"
  echo "  --image                image name (default: regenie)"
  echo "  --push                 push the built image after a successful build (default: build only, no push)"
  echo "  --registry             refinery or sandbox (default: sandbox)"
  echo "  --base-regenie-docker  use this pre-built image instead of building regenie from scratch (default: blank, i.e. build from scratch)"
  exit 1
}

IMAGE="regenie"
FINNGEN_TAG=""
PUSH=false
REGISTRY="sandbox"
BASE_REGENIE_DOCKER=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --image) IMAGE="$2"; shift 2 ;;
    --version) FINNGEN_TAG="$2"; shift 2 ;;
    --push) PUSH=true; shift ;;
    --registry) REGISTRY="$2"; shift 2 ;;
    --base-regenie-docker) BASE_REGENIE_DOCKER="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1" >&2; usage ;;
  esac
done

[[ -z "$FINNGEN_TAG" ]] && { echo "--version is required" >&2; usage; }
if [[ "$REGISTRY" != "refinery" && "$REGISTRY" != "sandbox" ]]; then
  echo "--registry must be 'refinery' or 'sandbox'" >&2
  usage
fi

# Registry path constants
REFINERY_PATH="europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/"
SB_PATH="eu.gcr.io/finngen-sandbox-v3-containers/"

BASE_PATH="$REFINERY_PATH"
[[ "$REGISTRY" == "sandbox" ]] && BASE_PATH="$SB_PATH"
REPO="${BASE_PATH}${IMAGE}"

#change directory to parent of this script so docker file can always be found without having to worry from which path you launch the script
APP_ROOT="$(dirname "$(dirname "$(readlink -fm "$0")")")"
cd $APP_ROOT

## which htslib version to use for tabix etc.
HTSLIB_VER=1.14

# BASE_REGENIE_DOCKER is set from --base-regenie-docker above; blank (the default) builds regenie from scratch

# only re-clone if local master doesn't match upstream master -- avoids re-cloning on every run while
# still catching the actual bug that started this: git silently no-op'ing on an existing directory
# (fatal: destination path already exists) left a 2+ year stale source tree with nothing catching it.
# Checked unconditionally (even when --base-regenie-docker skips the from-scratch build below) since
# TAG is always derived from regenie's own VERSION file.
# Note: this directory lives under the Dropbox-synced repo path. If Dropbox's sync daemon is actively
# touching it (observed once as multi-minute hangs on plain git/rm operations here), a re-clone below
# may hang too -- if that happens, pause Dropbox sync (or set up selective sync excluding this path).
REGENIE_SRC_DIR="$APP_ROOT/regenie"
REGENIE_REPO_URL="https://github.com/rgcgithub/regenie.git"

REMOTE_HASH=$(git ls-remote "$REGENIE_REPO_URL" HEAD | cut -f1)
if [[ -z "$REMOTE_HASH" ]]; then
  echo "Could not reach $REGENIE_REPO_URL to check for updates -- not attempting the build" >&2
  exit 1
fi

LOCAL_HASH=""
[[ -d "$REGENIE_SRC_DIR/.git" ]] && LOCAL_HASH=$(git -C "$REGENIE_SRC_DIR" rev-parse HEAD 2>/dev/null)

if [[ "$LOCAL_HASH" == "$REMOTE_HASH" ]]; then
  echo "Local regenie clone already at latest master ($REMOTE_HASH), skipping re-clone"
else
  echo "Local regenie clone stale or missing (local=${LOCAL_HASH:-none}, remote=$REMOTE_HASH) -- re-cloning into $REGENIE_SRC_DIR ..."
  rm -rf "$REGENIE_SRC_DIR"
  git clone "$REGENIE_REPO_URL" "$REGENIE_SRC_DIR"
  if [[ $? -ne 0 ]]; then
    echo "git clone failed -- not attempting the build" >&2
    exit 1
  fi
fi
cd "$REGENIE_SRC_DIR"

#optional
## uses master by default but comment out to use specific version
#git checkout c1daf24

TAG=$(cat VERSION)"_"$FINNGEN_TAG

if [[ -z "$BASE_REGENIE_DOCKER" ]]; then
  IMG_NAME=regenie:v$(cat VERSION)".gz"
  echo "Building base regenie version " $(cat VERSION)
  make docker-build MKLROOT=1 STATIC=1 HAS_BOOST_IOSTREAM=1
  if [[ $? -ne 0 ]]; then
    echo "make docker-build failed -- not attempting the finngen layer build" >&2
    exit 1
  fi

  # give the bare regenie image its own registry-qualified tag (same tag portion as IMG_NAME, e.g.
  # v4.1.2.gz), separate from the finngen image's tag below -- so it can be pushed/reused directly
  # (e.g. via --base-regenie-docker) without needing to hand-tag it after the fact
  BARE_TAG="$REPO:${IMG_NAME#*:}"
  docker tag $IMG_NAME $BARE_TAG
  echo "Pushing bare regenie image to $BARE_TAG"
  docker push $BARE_TAG
else
  echo "Skipping building base regenie and using $BASE_REGENIE_DOCKER"
  IMG_NAME=$BASE_REGENIE_DOCKER
fi

cd "$APP_ROOT"

echo $PWD

# tag the finngen layer's output directly as the final repo:tag -- never reuse $IMG_NAME for this, or
# the bare regenie image it points to gets overwritten, and --base-regenie-docker can never cache-hit
# against it again (Docker's cache is keyed on that image having been a FROM before; a name that used
# to mean "bare regenie" now meaning "full finngen image" has no such history and forces a full rebuild)
echo "Building finngen regenie-pipeline docker based on $IMG_NAME and using htslib $HTSLIB_VER"
docker build -f docker/Dockerfile -t $REPO:$TAG --build-arg base_image=$IMG_NAME --build-arg HTSLIB_VER=$HTSLIB_VER  .
if [[ $? -ne 0 ]]; then
  echo "docker build failed -- not tagging or pushing" >&2
  exit 1
fi
echo "Tagged $REPO:$TAG"

if [[ "$PUSH" == true ]]; then
  echo "Pushing to docker repo $REPO:$TAG"
  docker push $REPO:$TAG
else
  echo "Skipping push (pass --push to also push)"
fi

echo
echo "=== Summary ==="
echo "Bare regenie image: ${BARE_TAG:-$BASE_REGENIE_DOCKER (reused, not built this run)}"
echo "Finngen image:      $REPO:$TAG"
