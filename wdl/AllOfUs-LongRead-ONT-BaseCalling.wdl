version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"

task fast5ToPod5 {
    input {
        File inputFAST5
        String pod5Docker = "erictdawson/pod5tools"
    }

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 4,
        "gbRAM": 16,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputFAST5, "GB")) + 80 else runtime_attributes.diskGB

    String outbase = basename(inputFAST5, ".fast5")
    command <<<
        pod5 convert fast5 --input ~{inputFAST5} ~{outbase}.pod5
    >>>

    output {
        File outputPOD5 = "~{outbase}.pod5"
    }

    runtime {
        docker : "~{pod5Docker}"
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

task Dorado {
    input {
        File inputPOD5
        File refTarball
        String model = "dna_r10.4.1_e8.2_400bps_hac@v4.1.0"

        String doradoDocker = "ontresearch/dorado:sha31bb7b1fa487fb5b78d765406ecb0aa5ab78ef09"
    }
    Int sort_threads = 4

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 24,
        "gbRAM": 120,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    GPUAttributes gpu_attributes = {
        "gpuModel" : "nvidia-tesla-t4",
        "nGPU" : 4,
        "gpuDriverVersion": "535.104.05"
    }
    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    # Int mm2_threads = if nThreads - sort_threads >= mapThreads then mapThreads else nThreads - sort_threads -1
    String outbase = basename(basename(basename(inputPOD5, ".gz"), ".fq"), ".fastq")
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputPOD5, "GB") * 3.2) + ceil(size(refTarball, "GB") * 3) + 80 else runtime_attributes.diskGB
    String localTarball = basename(refTarball)
    String ref = basename(refTarball, ".tar")

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        mv ~{refTarball} ~{localTarball} && \
        tar xvf ~{localTarball} && \
        dorado download --model ~{model} && \
        dorado basecaller ~{"--reference " + ref} ~{model}
    >>>
    output {
        File? outputFASTQ = "~{outbase}.fastq.gz"
        File? outputBAM = "~{outbase}.bam"
        File? outputBAI= "~{outbase}.bam.bai"
    }
    runtime {
        docker : "~{doradoDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : runtime_attributes.nThreads
        memory : "~{runtime_attributes.gbRAM} GB"
        hpcMemory : runtime_attributes.gbRAM
        hpcQueue : "~{runtime_attributes.hpcQueue}"
        hpcRuntimeMinutes : runtime_attributes.runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : runtime_attributes.maxPreemptAttempts

        gpuType : "~{gpu_attributes.gpuModel}"
        gpuCount : gpu_attributes.nGPU
        nvidiaDriverVersion : "~{gpu_attributes.gpuDriverVersion}"
    }
}


# workflow {
#     input {
#         File inputSquigglefile
#         File inputRefTarball
#     }

#     output {
        
#     }



# }