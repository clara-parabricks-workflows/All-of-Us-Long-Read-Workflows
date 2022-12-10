All of Us ONT Long Read Pipelines
---------------------

# Overview
This repository contains workflows for alignment of ONT long reads and variant
calling of SNVs, indels, SVs, and STRS in the [Workflow Description Language](https://github.com/openwdl/wdl).

# Usage
These workflows may be used locally, with Docker, or on Terra.

## Building Dockerfiles
The docker images used in this repository can be built like so:

```bash
cd docker

docker build -f <app>.Dockerfile -t <user>/<app> .

## For a concrete example:
docker build -f clair3.Dockerfile -t erictdawson/clair3
```
# Workflows

## Alignment
No alignment workflow has been implemented. We recommend [minimap2](https://github.com/lh3/minimap2) for alignment.

## Variant calling
The following tools are used for variant calling (filled boxes represented implemented callers):

- [X] clair3 (SNVs and indels)
- [X] DeepVariant (SNVs and indels) (with Parabricks)
- [X] sniffles2 (SVs)
- [ ] STRspy (STRs)