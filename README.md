All of Us ONT Long Read Pipelines
---------------------

# Overview
This repository contains workflows for alignment of ONT long reads and variant
calling of SNVs, indels, SVs, and STRS in the [Workflow Description Language](https://github.com/openwdl/wdl).

# Usage
These workflows may be used locally, with Docker, or on Terra.

Running locally with or without Docker can be done via Cromwell using the config files in the `configs` directory.

To download Cromwell:

```bash
## Update the version as needed
export version=81
wget https://github.com/broadinstitute/cromwell/releases/download/${version}/cromwell-${version}.jar
https://github.com/broadinstitute/cromwell/releases/download/${version}/womtool-${version}.jar
```

To run using Docker:

```bash
java -Dconfig.file=config/localDocker.wdl.conf -jar cromwell-81.jar run -i inputs.local.json wdl/AllOfUs-LongRead-ONT-AlignmentAndPhasing.wdl
```

To run using GCP, you can modify the gcp_template.wdl to include your project and billing ID / bucket and run:

```bash
java -Dconfig.file=config/gcp_template.wdl.conf -jar cromwell-81.jar run -i inputs.local.json wdl/AllOfUs-LongRead-ONT-AlignmentAndPhasing.wdl
```

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
The `wdl/AllOfUs-LongRead-ONT-AlignmentAndPhasing.wdl ` workflow runs minimap2 for alignment, sorting and indexing and then runs BAM tagging
with WhatsHap. The final output is a phased, tagged BAM file.

## Variant calling
The following tools are used for variant calling (filled boxes represented implemented callers):

- [X] clair3 (SNVs and indels)
- [X] DeepVariant (SNVs and indels) (with Parabricks)
- [X] sniffles2 (SVs)
- [ ] STRspy (STRs)