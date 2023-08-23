version 1.0

task fast5ToPod5 {
    input {
        File inputFAST5
    }
    
    String outbase = basename(inputFAST5, ".fast5")
    command <<<
        pod5 convert fast5 --input ~{inputFAST5} ~{outbase}.pod5
    >>>

    output {
        File outputPOD5 = "~{outbase}.pod5"
    }

    runtime {
        docker : "~{minimapDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

task Dorado {
    input {
        File inputPOD5
        File? refTarball

        Int nThreads = 24
        String doradoDocker = "ontresearch/dorado:sha31bb7b1fa487fb5b78d765406ecb0aa5ab78ef09"
        Int diskGB = 0
        Int nGPU = 2
        String gpuModel = "nvidia-tesla-t4"
        Int maxPreemptAttempts = 3
        Int runtimeMinutes = 300
        String hpcQueue = "norm"
        Int gbRAM = 87
    }
    Int sort_threads = 4
    ## Put a ceiling on mm2_threads so as not to oversubscribe our VM
    ## mm2_threads = min(mapThreads, nThreads - sort_threads - 1)
    Int mm2_threads = if nThreads - sort_threads >= mapThreads then mapThreads else nThreads - sort_threads -1
    String outbase = basename(basename(basename(inputFASTQ, ".gz"), ".fq"), ".fastq")
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputFASTQ, "GB") * 3.2) + ceil(size(inputReference, "GB") * 3) + 80 else diskGB
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
         -@ ~{sort_threads} - \
         > ~{outbase}.bam && \
         samtools index ~{outbase}.bam
    }
    output {
        File? outputFASTQ
        File? outputBAM = "~{outbase}.bam"
        File? outputBAI= "~{outbase}.bam.bai"
    }
    runtime {
        docker : "~{doradoDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        gpuType : "~{gpuModel}"
        gpuCount : nGPU
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}


workflow {
    input {
        File inputSquigglefile
        File inputRefTarball
    }

    output {
        
    }



}