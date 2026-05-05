#!/usr/bin/env bash
# ============================================================
#  build_and_push_images.sh
#  Builds both custom Docker images and pushes them to Docker Hub.
#
#  Usage:
#    export DOCKERHUB_USER=youruser
#    export DOCKERHUB_PASS=yourpassword   # or use 'docker login' beforehand
#    bash build_and_push_images.sh
# ============================================================
set -euo pipefail

DOCKERHUB_USER="${DOCKERHUB_USER:-youruser}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DOCKER_DIR="${SCRIPT_DIR}/docker"

# ---------------------------------------------------------------
# Image 1: bwa-samtools
# ---------------------------------------------------------------
BWA_TAG="${DOCKERHUB_USER}/bwa-samtools:1.0.0"
echo ">>> Building ${BWA_TAG} ..."
docker build \
    --platform linux/amd64 \
    -t "${BWA_TAG}" \
    -f "${DOCKER_DIR}/Dockerfile.bwa-samtools" \
    "${DOCKER_DIR}"

docker push "${BWA_TAG}"

# ---------------------------------------------------------------
# Image 2: r-ggplot2
# ---------------------------------------------------------------
R_TAG="${DOCKERHUB_USER}/r-ggplot2:1.0.0"
echo ""
echo ">>> Building ${R_TAG} ..."
docker build \
    --platform linux/amd64 \
    -t "${R_TAG}" \
    -f "${DOCKER_DIR}/Dockerfile.r-ggplot2" \
    "${DOCKER_DIR}"

echo ">>> Pushing ${R_TAG} ..."
docker push "${R_TAG}"

echo ""
echo "=== All images built and pushed successfully ==="
echo "  ${BWA_TAG}"
echo "  ${R_TAG}"
echo ""
echo "Update nextflow.config container profile with your Docker Hub username:"
echo "  sed -i 's/youruser/${DOCKERHUB_USER}/g' nextflow.config"
