version 1.0

task AlignBam {
    input {
        File inputFASTQ
        File inputReference
        String sampleName
        String mm2Preset
        Int nThreads = 32
        Int mapThreads = 28
    }
    Int sort_threads = 4
    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    Int mm2_threads = if nThreads - sort_threads >= mapThreads then mapThreads else nThreads - sort_threads -1
    
    String outbase = basename(basename(basename(inputFASTQ, ".gz"), ".fq"), ".fastq")
    command {
        minimap2 \
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
         -@ ~{sort_threads} - \
         > ~{outbase}.bam && \
         samtools index ~{outbase}.bam

    }
    output {
        File outputBAM = "~{outbase}.bam"
        File outputBAI= "~{outbase}.bam.bai"
    }
    runtime {
        docker : "~{minimapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        # gpuType : "~{gpuModel}"
        # gpuCount : nGPU
        # nvidiaDriverVersion : "~{gpuDriverVersion}"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

task phaseAndTagBam {
    input {
        File inputBAM
        File inputBAI
        File inputVCF
        File inputRefTarball

        String whatshapDocker = "erictdawson/whatshap"
    }

    String outbase = basename(inputBAM, ".bam")
    String localRefTarball = basename(inputBAM)
    String inputReference = basename(inputRefTarball, ".tar")
    command {
        mv ~{inputRefTarball} ~{localRefTarball} && \
        tar xvf ~{localRefTarball} && \
        whatshap haplotag \
        -o ~{outbase}.phased.bam \
        -r ~{inputReference} \
        PHASED_SNPS.vcf.gz ~{inputBAM} && \
        samtools index ~{outbase}.phased.bam
    }
    output {
        File outputTaggedBam = "~{outbase}.phased.bam"
        File outputTaggedBAI = "~{outbase}.phased.bam.bai"
    }
    runtime {
        docker : "~{whatshapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        # gpuType : "~{gpuModel}"
        # gpuCount : nGPU
        # nvidiaDriverVersion : "~{gpuDriverVersion}"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

workflow AoU_ONT_Alignment {
    input {
        File inputFASTQ
        File reference
    }

    call minimap2 {
        input:
    }

    call phaseAndTagBam {
        input:
    }

    output {

    }
}