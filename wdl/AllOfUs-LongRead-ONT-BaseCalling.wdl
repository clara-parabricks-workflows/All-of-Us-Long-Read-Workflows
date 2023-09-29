version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"

task fast5ToPod5 {
    input {
        File inputFAST5tarball
        String pod5Docker = "erictdawson/pod5tools"
    }

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 8,
        "gbRAM": 26,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputFAST5tarball, "GB") * 4) + 80 else runtime_attributes.diskGB

    String outbase = basename(inputFAST5tarball, ".tar")
    String outputDir = outbase + ".pod5"
    String outputTarball = outputDir + ".tar"
    command <<<
        mkdir ~{outbase}.pod5s && \
        tar xf ~{inputFAST5tarball} && \
        pod5 convert fast5 --one-to-one ~{outbase}/ --input ~{outbase}/ ~{outputDir} && \
        tar cf ~{outputTarball} ~{outputDir}
    >>>

    output {
        File outputPOD5tarball = "~{outputTarball}"
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

# task splitPOD5ByChannel {
#     input {
#         File inputPOD5tarball
#     }

#     RuntimeAttributes runtime_attributes = {
#         "diskGB": 0,
#         "nThreads": 24,
#         "gbRAM": 120,
#         "hpcQueue": "norm",
#         "runtimeMinutes": 600,
#         "maxPreemptAttempts": 3,
#     }

#     String outbase = basename(inputPOD5tarball, ".tar")

#     command <<<
#     tar xf ~{inputPOD5tarball} && \
#     time pod5 view ~{outbase} --include "read_id,channel" --output summary.tsv && \
#     time pod5 subset ~{pod5dir} --summary summary.tsv --columns channel --output split_by_channel
#     >>>

#     output {

#     }

#     runtime {

#     }
# }

# task DoradoNoAlign {
#     input {

#     }
#     command <<<
#     >>>
#     output {

#     }
#     runtime {

#     }
# }

task DoradoWithAlignment {
    input {
        File inputPOD5tarball
        File inputRefTarball
        String model = "dna_r10.4.1_e8.2_400bps_hac@v4.1.0"

        String doradoDocker = "ontresearch/dorado"
    }

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

    Int sort_threads = 6

    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    # Int mm2_threads = if nThreads - sort_threads >= mapThreads then mapThreads else nThreads - sort_threads -1
    String outbase = basename(inputPOD5tarball, ".tar")
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputPOD5tarball, "GB") * 3.2) + ceil(size(inputRefTarball, "GB") * 3) + 80 else runtime_attributes.diskGB
    String localTarball = basename(inputRefTarball)
    String ref = basename(inputRefTarball, ".tar")

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        mv ~{inputRefTarball} ~{localTarball} && \
        tar xvf ~{localTarball} && \
        tar xf ~{inputPOD5tarball} && \
        dorado download --model ~{model} && \
        dorado basecaller ~{"--reference " + ref} ~{model} ~{outbase} | \
        samtools sort --threads ~{sort_threads} -m12g -O BAM -o ~{outbase}.bam --write-index && \
        samtools index ~{outbase}.bam
    >>>

    output {
        File outputBAM = "~{outbase}.bam"
        File outputBAI= "~{outbase}.bam.bai"
        File outputCSI = "~{outbase}.bam.csi"
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


workflow DoradoBasecall {
    input {
        File inputFAST5tarball
        File inputRefTarball
        String pod5Docker = "erictdawson/pod5tools"
    }

    call fast5ToPod5{
        input:
            inputFAST5tarball=inputFAST5tarball,
            pod5Docker=pod5Docker
    }

    call DoradoWithAlignment{
        input:
            inputPOD5tarball=fast5ToPod5.outputPOD5tarball,
            inputRefTarball=inputRefTarball
    }


    output {
        File? outputPOD5tarball = fast5ToPod5.outputPOD5tarball
        File? outputBAM = DoradoWithAlignment.outputBAM
        File? outputBAI = DoradoWithAlignment.outputBAI        
    }
}