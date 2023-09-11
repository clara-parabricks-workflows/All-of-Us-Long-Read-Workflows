version 1.0

task phaseVCF {
    input {
        File inputBAM
        File inputBAI
        File inputVCF
        File inputTBI
        File inputRefTarball

        String whatshapDocker = "erictdawson/whatshap"
        Int diskGB = 0
        Int nThreads = 24
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        Int maxPreemptAttempts = 3
    }

    String outbase = basename(basename(inputVCF, ".gz"), ".vcf")
    String localRefTarball = basename(inputBAM)
    String inputReference = basename(inputRefTarball, ".tar")

    command {
        whatshap phase \
        -o ~{outbase}.phased.vcf \
        --reference ~{inputReference} \
        ~{inputVCF} \
        ~{inputBAM} && \
        bgzip ~{outbase}.phased.vcf && \
        tabix ~{outbase}.phased.vcf.gz
    }
    output {
        File outputPhasedVCF = "~{outbase}.phased.vcf.gz"
        File outputPhasedTBI = "~{outbase}.phased.vcf.gz.tbi"
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

task haplotagBAM {
    input {
        File inputBAM
        File inputBAI
        File inputVCF
        File inputTBI
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
        --output-threads ~{nThreads}
        -o ~{outbase}.phased.bam \
        -r ~{inputReference} \
        ~{inputVCF} \
        ~{inputBAM} && \
        samtools index ~{outbase}.phased.bam
    }

    output {
        File outputHaplotaggedBAM = "~{outbase}.phased.bam"
        File outputHaplotaggedBAI = "~{outbase}.phased.bam.bai"
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

workflow AoU_ONT_Phase {
    input {
        File inputBAM
        File inputBAI
        File inputRefTarball
        File inputVCF
        File inputTBI
    }

    call phaseVCF {
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputVCF=inputVCF,
            inputTBI=inputTBI
    }

    call haplotagBAM {
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputVCF=phaseVCF.outputPhasedVCF,
            inputTBI=phaseVCF.outputPhasedTBI
    }

    output {
        File phasedVCF = phaseVCF.outputPhasedVCF
        File phasedTBI = phaseVCF.outputPhasedTBI
        File haplotaggedBAM = haplotagBAM.outputHaplotaggedBAM
        File haplotaggedBAI = haplotagBAM.outputHaplotaggedBAI
    }

}