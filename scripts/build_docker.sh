#!/bin/bash
set -u

usage() {
  echo "Usage: $0 --version <finngen_tag> [--image <name>] [--push] [--registry refinery|sandbox] [--base-regenie-docker <image>]"
  echo "  --version              required. Suffix appended to regenie's own VERSION file to build the pushed tag (e.g. cond_firth)"
  echo "  --image                image name (default: regenie)"
  echo "  --push                 also push the finngen image after a successful build (the bare regenie image is always pushed)"
  echo "  --registry             refinery or sandbox (default: sandbox)"
  echo "  --base-regenie-docker  use this pre-built image instead of building regenie from scratch"
  exit 1
}

die() { echo "$*" >&2; exit 1; }
run() { "$@" || die "failed: $*"; }

IMAGE="regenie"; FINNGEN_TAG=""; PUSH=false; REGISTRY="sandbox"; BASE_REGENIE_DOCKER=""

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
[[ "$REGISTRY" == "refinery" || "$REGISTRY" == "sandbox" ]] || { echo "--registry must be 'refinery' or 'sandbox'" >&2; usage; }

REFINERY_PATH="europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/"
SB_PATH="eu.gcr.io/finngen-sandbox-v3-containers/"
[[ "$REGISTRY" == "sandbox" ]] && REPO="${SB_PATH}${IMAGE}" || REPO="${REFINERY_PATH}${IMAGE}"

# repo root, so the build context is found regardless of launch dir
APP_ROOT="$(dirname "$(dirname "$(readlink -fm "$0")")")"
cd "$APP_ROOT"

HTSLIB_VER=1.14
REGENIE_SRC_DIR="$APP_ROOT/regenie"
REGENIE_REPO_URL="https://github.com/rgcgithub/regenie.git"

# re-clone only if stale (a plain `git clone` into an existing dir silently no-ops otherwise)
REMOTE_HASH=$(git ls-remote "$REGENIE_REPO_URL" HEAD | cut -f1)
[[ -n "$REMOTE_HASH" ]] || die "could not reach $REGENIE_REPO_URL to check for updates"
LOCAL_HASH=""
[[ -d "$REGENIE_SRC_DIR/.git" ]] && LOCAL_HASH=$(git -C "$REGENIE_SRC_DIR" rev-parse HEAD 2>/dev/null)
if [[ "$LOCAL_HASH" != "$REMOTE_HASH" ]]; then
  echo "regenie clone stale or missing (local=${LOCAL_HASH:-none}, remote=$REMOTE_HASH) -- re-cloning"
  rm -rf "$REGENIE_SRC_DIR"
  run git clone "$REGENIE_REPO_URL" "$REGENIE_SRC_DIR"
else
  echo "regenie clone already at latest master ($REMOTE_HASH)"
fi
cd "$REGENIE_SRC_DIR"

TAG="$(cat VERSION)_${FINNGEN_TAG}"

# if the expected bare-regenie tag already exists remotely, offer to skip rebuilding it entirely
if [[ -z "$BASE_REGENIE_DOCKER" ]]; then
  EXPECTED_TAG="$REPO:v$(cat VERSION).gz"
  if docker manifest inspect "$EXPECTED_TAG" >/dev/null 2>&1; then
    read -r -p "$EXPECTED_TAG already exists remotely -- rebuild it anyway? [y/N] " REPLY
    [[ "$REPLY" =~ ^[Yy]$ ]] || BASE_REGENIE_DOCKER="$EXPECTED_TAG"
  fi
fi

if [[ -z "$BASE_REGENIE_DOCKER" ]]; then
  IMG_NAME="regenie:v$(cat VERSION).gz"
  echo "Building base regenie $(cat VERSION)"
  run make docker-build MKLROOT=1 STATIC=1 HAS_BOOST_IOSTREAM=1

  # own tag, pushed immediately, so --base-regenie-docker can reuse it and the finngen build below never overwrites it
  BARE_TAG="$REPO:${IMG_NAME#*:}"
  run docker tag "$IMG_NAME" "$BARE_TAG"
  run docker push "$BARE_TAG"
else
  echo "Using existing base image $BASE_REGENIE_DOCKER"
  IMG_NAME="$BASE_REGENIE_DOCKER"
fi

cd "$APP_ROOT"
echo "Building finngen layer on $IMG_NAME (htslib $HTSLIB_VER)"
run docker build -f docker/Dockerfile -t "$REPO:$TAG" --build-arg base_image="$IMG_NAME" --build-arg HTSLIB_VER="$HTSLIB_VER" .

if [[ "$PUSH" == true ]]; then
  run docker push "$REPO:$TAG"
fi

echo
echo "=== Summary ==="
echo "Bare regenie image: ${BARE_TAG:-$BASE_REGENIE_DOCKER (reused, not built this run)}"
if [[ "$PUSH" == true ]]; then
  echo "Finngen image:      $REPO:$TAG"
else
  echo "Finngen image:      $REPO:$TAG (not pushed -- pass --push)"
fi
