version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/vcf.wdl" as vcf
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/sam.wdl" as sam

task phaseVCF {
    input {
        File inputBAM
        File inputBAI
        File inputVCF
        File inputTBI
        File inputRefTarball
        String? sampleName

        String whatshapDocker = "erictdawson/whatshap"
        Int diskGB = 0
        Int nThreads = 24
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        Int maxPreemptAttempts = 3
    }

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 8,
        "gbRAM": 59,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    String outbase = basename(basename(inputVCF, ".gz"), ".vcf")
    String localRefTarball = basename(inputBAM)
    String ref = basename(inputRefTarball, ".tar")
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(inputRefTarball, "GB")) + ceil(size(inputBAI, "GB")) + 65 else runtime_attributes.diskGB

    command {
        cp ~{inputRefTarball} ~{localRefTarball} && \
        tar vxf ~{localRefTarball} && \
        whatshap phase \
        -o ~{outbase}.phased.vcf \
        ~{"--sample " + sampleName} \
        --reference ~{ref} \
        ~{inputVCF} \
        ~{inputBAM}
    }
    output {
        File outputPhasedVCF = "~{outbase}.phased.vcf"
    }
    runtime {
        docker : "~{whatshapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{runtime_attributes.gbRAM} GB"
        hpcMemory : runtime_attributes.gbRAM
        hpcQueue : "~{runtime_attributes.hpcQueue}"
        hpcRuntimeMinutes : runtime_attributes.runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : runtime_attributes.maxPreemptAttempts
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
    }

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 8,
        "gbRAM": 59,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    String outbase = basename(inputBAM, ".bam")
    String localRefTarball = basename(inputRefTarball)
    String ref= basename(localRefTarball, ".tar")
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(inputRefTarball, "GB") * 3) + 80 else runtime_attributes.diskGB
    
    command {
        cp ~{inputRefTarball} ~{localRefTarball} && \
        tar xvf ~{localRefTarball} && \
        whatshap haplotag \
        --output-threads ~{runtime_attributes.nThreads} \
        --output ~{outbase}.phased.bam \
        --reference ~{ref} \
        ~{inputVCF} \
        ~{inputBAM}
    }

    output {
        File outputHaplotaggedBAM = "~{outbase}.phased.bam"
    }
    runtime {
        docker : "~{whatshapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : runtime_attributes.nThreads
        memory : "~{runtime_attributes.gbRAM} GB"
        hpcMemory : runtime_attributes.gbRAM
        hpcQueue : "~{runtime_attributes.hpcQueue}"
        hpcRuntimeMinutes : runtime_attributes.runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : runtime_attributes.maxPreemptAttempts
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
            inputRefTarball=inputRefTarball,
            inputBAI=inputBAI,
            inputVCF=inputVCF,
            inputTBI=inputTBI
    }

    RuntimeAttributes compress_attributes = {
        "diskGB": 0,
        "nThreads": 4,
        "gbRAM": 9,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    call vcf.compressAndIndexVCF as compress_phased_vcf {
        input:
            inputVCF=phaseVCF.outputPhasedVCF,
            attributes=compress_attributes
    }

    call haplotagBAM {
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputRefTarball=inputRefTarball,
            inputVCF=compress_phased_vcf.outputVCFGZ,
            inputTBI=compress_phased_vcf.outputTBI
    }

    call sam.indexBAM as indexTaggedBAM {
        input:
            inputBAM=haplotagBAM.outputHaplotaggedBAM
    }



    output {
        File phasedVCF = compress_phased_vcf.outputVCFGZ
        File phasedTBI = compress_phased_vcf.outputTBI
        File haplotaggedBAM = haplotagBAM.outputHaplotaggedBAM
        File haplotaggedBAI = indexTaggedBAM.outputBAI
    }

}