#!/bin/bash
set -e

VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo "dev")
IMAGE_NAME="fusilli_${VERSION}.sif"

singularity build "${IMAGE_NAME}" fusilli.def
echo "Built: ${IMAGE_NAME}"
