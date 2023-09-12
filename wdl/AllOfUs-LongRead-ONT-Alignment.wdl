version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"

task AlignBam {
    input {
        File inputFASTQ
        File inputReference
        File? referenceIndex
        String sampleName
        String mm2Preset = "map-ont"
        Int mapThreads = 28
        Int sortThreads = 4
        String minimapDocker = "erictdawson/minimap2"

        RuntimeAttributes attributes = {
            "diskGB": 0,
            "nThreads": 64,
            "gbRAM": 126,
            "hpcQueue": "norm",
            "runtimeMinutes": 600,
            "gpuDriverVersion": "535.104.05",
            "maxPreemptAttempts": 3,
            "zones": ["us-central1-a", "us-central1-b", "us-central1-c"]
        }
    }
    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    Int mm2_threads = if attributes.nThreads - sortThreads >= mapThreads then mapThreads else attributes.nThreads - sortThreads -1
    String outbase = basename(basename(basename(inputFASTQ, ".gz"), ".fq"), ".fastq")
    Int auto_diskGB = if attributes.diskGB == 0 then ceil(size(inputFASTQ, "GB") * 3.2) + ceil(size(inputReference, "GB") * 3) + 80 else attributes.diskGB
    command {
        time minimap2 \
        -Y \
        -H \
        -y \
        --MD \
        -t ~{mm2_threads} \
        -R "@RG\tSM:~{sampleName}\tID:~{sampleName}" \
        -ax ~{mm2Preset} \
        ~{inputReference} \
        ~{inputFASTQ} | \
        samtools sort \
         -m 6G \
         -@ ~{sortThreads} - \
         > ~{outbase}.bam && \
         samtools index -@4 ~{outbase}.bam

    }
    output {
        File outputBAM = "~{outbase}.bam"
        File outputBAI= "~{outbase}.bam.bai"
    }
    runtime {
        docker : "~{minimapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : attributes.nThreads
        memory : "~{attributes.gbRAM} GB"
        hpcMemory : attributes.gbRAM
        hpcQueue : "~{attributes.hpcQueue}"
        hpcRuntimeMinutes : attributes.runtimeMinutes
        zones : attributes.zones
        preemptible : attributes.maxPreemptAttempts
    }
}

workflow AoU_ONT_Alignment {
    input {
        File inputFASTQ
        File inputReference
        File? referenceIndex
        String sampleName = "sample"
        String mm2Preset = "map-ont"
        Int mapThreads = 32
        String minimapDocker = "erictdawson/minimap2"
        Int minimap_RAM = 62
        String minimap_queue = "norm"
        Int minimap_runtime_max = 600
    }

        RuntimeAttributes attributes = {
        "diskGB": 0,
        "nThreads": 64,
        "gbRAM": 125,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "gpuDriverVersion": "535.104.05",
        "maxPreemptAttempts": 3,
        "zones": ["us-central1-a", "us-central1-b", "us-central1-c"]
    }

    call AlignBam {
        input:
            inputFASTQ=inputFASTQ,
            inputReference=inputReference,
            referenceIndex=referenceIndex,
            sampleName=sampleName,
            mm2Preset=mm2Preset,
            minimapDocker=minimapDocker,
            attributes=attributes
    }


    output {
        File alignedBAM = AlignBam.outputBAM
        File alignedBAI = AlignBam.outputBAI
    }
}
