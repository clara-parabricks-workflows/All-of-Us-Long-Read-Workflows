All of Us Long Read Pipelines
---------------------

# Overview
This repository contains workflows for alignment of ONT long reads and variant
calling of SNVs, indels, SVs, and STRS in the [Workflow Description Language](https://github.com/openwdl/wdl).

# Usage

These workflows are available on [Dockstore]() and [Terra]().

## Running with existing POD5 / FAST5 / FASTQ / BAM / CRAM files on the cloud

See the [Whole Genome HG002 Case Study](docs/Whole-genome-HG002-case-study.md).

## Running locally and on other clouds

These workflows may be used locally, with Docker, or on Terra. They can also be modified to run on other clouds or HPC systems
with minimal effort.

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

## Basecalling

- [X] Pod5 Convert FAST5
- [X] Dorado (with alignment)

The basecalling workflow is available at the following links: [Dockstore](https://dockstore.org/my-workflows/github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/AllOfUs-LongRead-ONT-BaseCalling) [Github](https://github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/blob/main/wdl/AllOfUs-LongRead-ONT-BaseCalling.wdl) [Terra](https://app.terra.bio/#workspaces/clara-terra/All%20of%20US%20ONT)

## Alignment
The `wdl/AllOfUs-LongRead-ONT-Alignment.wdl ` workflow runs minimap2 for alignment, sorting and indexing. The final output is BAM file.

- [X] FASTQ -> BAM + BAI with minimap2 [Dockstore](https://dockstore.org/workflows/github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/AllOfUs-LongRead-ONT-Alignment:main?tab=info)

The alignment workflow can be found at the following links: [Dockstore](https://dockstore.org/my-workflows/github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/AllOfUs-LongRead-ONT-Alignment) [Github](https://github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/tree/main/wdl/AllOfUs-LongRead-ONT-Alignment.wdl) [Terra](https://app.terra.bio/#workspaces/clara-terra/All%20of%20US%20ONT)



## Phasing

- [X] Whatshap phaseVCF
- [X] Whatshap haplotag BAM

The phasing workflow can be found at the following links: [Dockstore](https://dockstore.org/my-workflows/github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/AllOfUs-LongRead-ONT-PhaseBam) [Github](https://github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/blob/main/wdl/AllOfUs-LongRead-ONT-PhaseBam.wdl) [Terra](https://app.terra.bio/#workspaces/clara-terra/All%20of%20US%20ONT)

## Variant calling
The following tools are used for variant calling (filled boxes represented implemented callers):

- [X] clair3 (SNVs and indels)
- [X] Parabricks DeepVariant (SNVs and indels)
- [X] sniffles2 (SVs)
- [X] spectre (CNVs)

The variant calling workflow can be found at the following links: [Dockstore](https://dockstore.org/my-workflows/github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/AllOfUs-LongRead-ONT-VariantCalling) [Github](https://github.com/clara-parabricks-workflows/All-of-Us-Long-Read-Workflows/blob/main/wdl/AllOfUs-LongRead-ONT-VariantCalling.wdl) [Terra](https://app.terra.bio/#workspaces/clara-terra/All%20of%20US%20ONT)