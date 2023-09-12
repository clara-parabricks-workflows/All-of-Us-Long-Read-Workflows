All of Us Long Read Pipelines
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
export version=84
wget https://github.com/broadinstitute/cromwell/releases/download/${version}/cromwell-${version}.jar
https://github.com/broadinstitute/cromwell/releases/download/${version}/womtool-${version}.jar
```

To validate a WDL:

```bash
make validate
```

To generate example inputs:

```bash
make inputs
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
The docker images used in this repository can be built like so, where user is a Dockerhub username and app corresponds to the basename of a dockerfile:

```bash
cd docker
make build APP=<app name> USER=<username>

## Here's an example:
cd docker
make build APP=spectre USER=eric
```

# Workflows

## Alignment
The `wdl/AllOfUs-LongRead-ONT-Alignment.wdl ` workflow runs minimap2 for alignment, sorting and indexing. The final output is BAM file.

- [X] FASTQ -> BAM + BAI with minimap2

## Phasing

- [ ] Phasing with WhatsHap


## Variant calling
The following tools are used for variant calling (filled boxes represented implemented callers):

- [X] clair3 (SNVs and indels)
- [X] DeepVariant (SNVs and indels) (with Parabricks)
- [X] sniffles2 (SVs)
- [ ] STRspy (STRs)
