#!/bin/bash
git clone git@github.com:rgcgithub/regenie.git

cd regenie

#optional
## uses master by default but omment out to use specific version
#git checkout c1daf24

## which htslib version to use for tabix etc.
HTSLIB_VER=1.14

## name to be appended to the resulting image name
FINNGEN_TAG=$1

### if this is given then does not rebuild regenie docker but uses this as base docker for building finngen docker
## leave empty or undefined to build regenie first.
BASE_REGENIE_DOCKER="eu.gcr.io/finngen-refinery-dev/regenie:v2.2.4.fixed.gz.mkl"

TAG=$(cat VERSION)"_"$FINNGEN_TAG

if [[ -z "$BASE_REGENIE_DOCKER" ]]; then
  IMG_NAME=regenie:v$(cat VERSION)".gz"
  echo "Building base regenie version " $(cat VERSION)
  make docker-build MKLROOT=1 STATIC=1 HAS_BOOST_IOSTREAM=1
else
  echo "Skipping building base regenie and using $BASE_REGENIE_DOCKER"
  IMG_NAME=$BASE_REGENIE_DOCKER
fi

cd ..

echo "Building finngen regenie-pipeline docker based on $IMG_NAME and using htslib $HTSLIB_VER"

docker build -f docker/Dockerfile -t $IMG_NAME --build-arg base_image=$IMG_NAME --build-arg HTSLIB_VER=$HTSLIB_VER .

echo "Pushing to docker repo eu.gcr.io/finngen-refinery-dev/regenie:$TAG"
docker tag $IMG_NAME eu.gcr.io/finngen-refinery-dev/regenie:$TAG
docker push eu.gcr.io/finngen-refinery-dev/regenie:$TAG
