version 1.0

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

workflow AoU_ONT_Phase {
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

}