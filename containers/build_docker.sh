#!/bin/bash
set -e

VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo "dev")
IMAGE_NAME="fusilli:${VERSION}"

docker build -t "${IMAGE_NAME}" -t "fusilli:latest" -f Dockerfile ..
echo "Built: ${IMAGE_NAME}"
