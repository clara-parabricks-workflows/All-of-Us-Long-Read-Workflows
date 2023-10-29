version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/sam.wdl"

task generatePOD5Summary {
    input {
        Array[File] inputPOD5s
        String pod5Docker = "erictdawson/pod5tools"
        RuntimeAttributes runtime_attributes = {
            "diskGB": 0,
            "nThreads": 8,
            "gbRAM": 26,
            "hpcQueue": "norm",
            "runtimeMinutes": 600,
            "maxPreemptAttempts": 3,
        }
    }
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputPOD5s, "GB") * 3.8) else runtime_attributes.diskGB

    command <<<
        pod5 view --threads ~{runtime_attributes.nThreads - 1} --include "read_id, channel" --output summary.tsv ~{sep=" " inputPOD5s}
    >>>
    output {
        File summary = "summary.tsv"
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

task fast5ToPod5 {
    input {
        File inputFAST5
        String pod5Docker = "erictdawson/pod5tools"
        
        RuntimeAttributes runtime_attributes = {
            "diskGB": 0,
            "nThreads": 2,
            "gbRAM": 9,
            "hpcQueue": "norm",
            "runtimeMinutes": 600,
            "maxPreemptAttempts": 3,
        }
    }

    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputFAST5, "GB") * 4) + 80 else runtime_attributes.diskGB

    String outbase = basename(inputFAST5, ".fast5")
    command <<<
        pod5 convert fast5 --output ~{outbase}.pod5  ~{inputFAST5}
    >>>

    output {
        File outputPOD5= "~{outbase}.pod5"
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

task splitPOD5ByChannel {
    input {
        Array[File] inputPOD5s
        File summaryFile
        String pod5Docker = "erictdawson/pod5tools"
    }

    RuntimeAttributes runtime_attributes = {
        "diskGB": 0,
        "nThreads": 32,
        "gbRAM": 84,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }


    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputPOD5s, "GB") * 3.8) else runtime_attributes.diskGB

    command <<<
        set -e
        mkdir split_by_channel && \
        time pod5 subset \
        --threads ~{runtime_attributes.nThreads - 2} \
        --summary ~{summaryFile} \
        --columns channel \
        --output split_by_channel \
        ~{sep=" " inputPOD5s}
    >>>

    output {
        Array[File] split_by_channel = glob("split_by_channel/*")
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
        String sampleName
        String? nameAnnotation
        Array[File] inputPOD5s
        File? inputRefTarball
        Boolean emitFASTQ = false
        Boolean duplex = false
        Int? minQScore
        Int? mm2_batch_size
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
    String outbase = if defined(nameAnnotation) then sampleName + "." + nameAnnotation else sampleName
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputPOD5s, "GB") * 4.8) + ceil(size(inputRefTarball, "GB") * 3) + 80 else runtime_attributes.diskGB
    String ref = if defined(inputRefTarball) then basename(select_first([inputRefTarball]), ".tar") else ""

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        ## Hack to create an input directory for FAST5 / POD5 files.
        mkdir -p data && \
        for i in `cat ~{write_lines(inputPOD5s)}`; do ln -s $i data/$(basename $i); done

        ~{if defined(inputRefTarball) then "tar xvf " + inputRefTarball + " -C `pwd` && " else ""} \
        dorado download --model ~{model} && \
        dorado \
            ~{if duplex then "duplex" else "basecaller"} \
            ~{"-I " + mm2_batch_size + "G"} \
            ~{if emitFASTQ then "--emit-fastq" else ""} \
            ~{"--min-qscore " + minQScore} \
            ~{if defined(inputRefTarball) then "--reference " + ref else ""} \
            ~{model} \
            data/ | \
        ~{if defined(inputRefTarball) then "samtools sort --threads " + sort_threads + " -m12g -O BAM -o " + outbase + ".bam" else " > " + outbase + ".fastq &&"} \
        ~{if defined(inputRefTarball) then "samtools index " + outbase + ".bam" else '"'}
    >>>

    output {
        File? outputBAM = "~{outbase}.bam"
        File? outputBAI= "~{outbase}.bam.bai"
        File? outputFASTQ = "~{outbase}.fastq"
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
        String sampleName
        String? nameAnnotation
        File inputList
        Boolean isPOD5 = true
        Boolean runSimplex = true
        Boolean runDuplex = false
        File? inputRefTarball
        Boolean emitFASTQ = false
        Int? minQScore
        String pod5Docker = "erictdawson/pod5tools"
    }

    Array[File] inputs = read_lines(inputList)

    if (!isPOD5){
        scatter(pod in inputs){
            call fast5ToPod5 {
                input:
                    inputFAST5=pod
            }
        }
    }

    Array[File] preparedPOD5s = select_first([fast5ToPod5.outputPOD5, inputs])

    if (runSimplex){
        call Dorado as simplex {
            input:
                sampleName=sampleName,
                nameAnnotation=nameAnnotation,
                inputPOD5s=preparedPOD5s,
                inputRefTarball=inputRefTarball,
                minQScore=minQScore,
                emitFASTQ=emitFASTQ,
                duplex=false
        }
    }

    if (runDuplex){
        call generatePOD5Summary{
            input:
                inputPOD5s=preparedPOD5s
        }
        call splitPOD5ByChannel {
            input:
                inputPOD5s=preparedPOD5s,
                summaryFile=generatePOD5Summary.summary
        }
        call Dorado as duplex {
            input:
                sampleName=sampleName,
                nameAnnotation=nameAnnotation,
                inputPOD5s=splitPOD5ByChannel.split_by_channel,
                inputRefTarball=inputRefTarball,
                minQScore=minQScore,
                emitFASTQ=emitFASTQ,
                duplex=true
        }
    }


    output {
        Array[File]? outputPOD5array = preparedPOD5s
        File? simplexBAM = simplex.outputBAM
        File? simplexBAI = simplex.outputBAI
        File? simplexFASTQ = simplex.outputFASTQ

        File? duplexBAM = duplex.outputBAM
        File? deplexBAI = duplex.outputBAI
        File? duplexFASTQ = duplex.outputFASTQ
    }
}