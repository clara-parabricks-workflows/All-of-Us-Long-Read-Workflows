version 1.0

task AlignBam {
    input {
        File inputFASTQ
        File inputRefTarball
        String sampleName
        String mm2Preset = "mapont"
        Int nThreads = 32
        Int mapThreads = 28

        String? minimapDocker = "erictdawson/minimap2"
        Int diskGB = 0
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        Int maxPreemptAttempts = 3
    }
    Int sort_threads = 4
    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    Int mm2_threads = if nThreads - sort_threads >= mapThreads then mapThreads else nThreads - sort_threads -1
    String outbase = basename(basename(basename(inputFASTQ, ".gz"), ".fq"), ".fastq")
    String inputReference = basename(inputRefTarball, ".tar")
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputFASTQ, "GB") * 3.2) + ceil(size(inputRefTarball, "GB") * 3) + 80 else diskGB
    command {
        tar xvf ~{inputRefTarball} && \
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
        Int diskGB = 0
        Int nThreads = 24
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        Int maxPreemptAttempts = 3
    }

    String outbase = basename(inputBAM, ".bam")
    String localRefTarball = basename(inputBAM)
    String inputReference = basename(inputRefTarball, ".tar")
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(inputRefTarball, "GB") * 3) + 80 else diskGB
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
        File inputRefTarball
        String sampleName = "sample"
        String mm2Preset = "mapont"
        Int mapThreads = 32
        String minimapDocker = "erictdawson/minimap2"
        Int minimap_RAM = 62
        String minimap_queue = "norm"
        Int minimap_runtime_max = 600

        File phaseVCF
    }

    call AlignBam {
        input:
            inputFASTQ=inputFASTQ,
            inputRefTarball=inputRefTarball,
            sampleName=sampleName,
            mm2Preset=mm2Preset,
            nThreads=mapThreads,
            minimapDocker=minimapDocker,
            gbRAM=minimap_RAM,
            hpcQueue=minimap_queue,
            runtimeMinutes=minimap_runtime_max,
            maxPreemptAttempts=3
    }

    call phaseAndTagBam {
        input:
            inputBAM=AlignBam.outputBAM,
            inputBAI=AlignBam.outputBAI,
            inputVCF=phaseVCF,
            inputRefTarball=inputRefTarball
    }

    output {
        File phasedAlignedBAM = phaseAndTagBam.outputTaggedBam
        File phasedAlignedBAI = phaseAndTagBam.outputTaggedBAI
    }
}