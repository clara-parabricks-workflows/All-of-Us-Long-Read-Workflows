version 1.0

task clair3 {
    input {
        File inputBAM
        File inputBAI
        File refTarball
        String platform = "ont"
        String modelPath = "/opt/models/ont"
        String sampleName = "Sample"
        File? targetsBed
        Boolean gvcfMode = false
        Boolean phaseMode = false
        
        Int nThreads = 24
        String clairDocker = "erictdawson/clair"
        Int diskGB = 0
        # Int nGPU = 2
        # String gpuModel = "nvidia-tesla-t4"
        # String gpuDriverVersion = "460.73.01"
        Int maxPreemptAttempts = 3
        Int runtimeMinutes = 300
        String hpcQueue = "norm"
        Int gbRAM = 87

    }
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(refTarball, "GB")) + ceil(size(inputBAI, "GB")) + 65 else diskGB
    String ref = basename(refTarball, ".tar")
    String outbase = basename(inputBAM, ".bam")
    command {
        /opt/bin/run_clair3.sh \
        --bam_fn=~{inputBAM} \   
        --ref_fn=~{ref} \      
        --threads=~{nThreads} \          
        --platform="ont" \            
        --model_path="~{modelPath}" \
        --output=~{outbase}.clair3  \
        ~{"--bed_fn=" + targetsBed} \
        ~{"--sample_name=" + sampleName} \
        ~{if gvcfMode then "--gvcf" else ""} \
        ~{if phaseMode then "--enable_phasing" else ""}
    }
    output {
        File pileupVCF = "~{outbase}.clair3/pileup.vcf.gz"
        File fullAlignmentVCF = "~{outbase}.clair3/full_alignment.vcf.gz"
        File mergeVCF = "~{outbase}.clair3/merge_output.vcf.gz"
    }
    runtime {
        docker : "~{clairDocker}"
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

task deepvariant {
    input {
        File inputBAM
        File inputBAI
        File inputRefTarball
        String pbPATH = "pbrun"
        File? pbLicenseBin
        String? pbDocker
        Boolean gvcfMode = false
        Int nGPU = 4
        String gpuModel = "nvidia-tesla-t4"
        String gpuDriverVersion = "460.73.01"
        Int nThreads = 24
        Int gbRAM = 120
        Int diskGB = 0
        Int runtimeMinutes = 600
        String hpcQueue = "gpu"
        Int maxPreemptAttempts = 3
    }
    String ref = basename(inputRefTarball, ".tar")
    String localTarball = basename(inputRefTarball)
    String outbase = basename(inputBAM, ".bam")

    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB")) + ceil(size(inputRefTarball, "GB")) + ceil(size(inputBAI, "GB")) + 65 else diskGB

    String outVCF = outbase + ".deepvariant" + (if gvcfMode then '.g' else '') + ".vcf"
    command {
        mv ~{inputRefTarball} ${localTarball} && \
        time tar xvf ~{localTarball} && \
        time ${pbPATH} deepvariant \
        ~{if gvcfMode then "--gvcf " else ""} \
        --ref ${ref} \
        --in-bam ${inputBAM} \
        --out-variants ~{outVCF} \
        ~{"--license-file " + pbLicenseBin} && \
        bgzip ~{outVCF} && \
        tabix ~{outVCF}.gz
    }
    output {
        File outputVCF = "~{outVCF}.gz"
        File outputTBI = "~{outVCF}.gz.tbi"
    }
    runtime {
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        gpuType : "~{gpuModel}"
        gpuCount : nGPU
        nvidiaDriverVersion : "~{gpuDriverVersion}"
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

# task mosdepth {
#     input {
# Int diskGB = 0
# Int nThreads = 24
# Int gbRAM = 62
# String hpcQueue = "norm"
# Int runtimeMinutes = 240
# String mosdepthDocker = ""
# Int maxPreemptAttempts = 3
#     }
#     command {
#     }
#     output {
#     }
#     runtime {
#         docker : "~{mosdepthDocker}"
#         disks : "local-disk ~{auto_diskGB} SSD"
#         cpu : nThreads
#         memory : "~{gbRAM} GB"
#         hpcMemory : gbRAM
#         hpcQueue : "~{hpcQueue}"
#         hpcRuntimeMinutes : runtimeMinutes
#         zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
#         preemptible : maxPreemptAttempts
#     }
# }

# task strspy {
#     input {
# Int diskGB = 0
# Int nThreads = 24
# Int gbRAM = 62
# String hpcQueue = "norm"
# Int runtimeMinutes = 240
# String strspyDocker = ""
# Int maxPreemptAttempts = 3
#     }
#     command {
#     }
#     output {
#     }
#     runtime {
#         docker : "~{strspyDocker}"
#         disks : "local-disk ~{auto_diskGB} SSD"
#         cpu : nThreads
#         memory : "~{gbRAM} GB"
#         hpcMemory : gbRAM
#         hpcQueue : "~{hpcQueue}"
#         hpcRuntimeMinutes : runtimeMinutes
#         zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
#         preemptible : maxPreemptAttempts
#     }
# }

# task phasing {
#     input {
# Int diskGB = 0
# Int nThreads = 24
# Int gbRAM = 62
# String hpcQueue = "norm"
# Int runtimeMinutes = 240
# String phasingDocker = ""
# Int maxPreemptAttempts = 3
#     }
#     command {
#     }
#     output {
#     }
#     runtime {
#         docker : "~{phasingDocker}"
#         disks : "local-disk ~{auto_diskGB} SSD"
#         cpu : nThreads
#         memory : "~{gbRAM} GB"
#         hpcMemory : gbRAM
#         hpcQueue : "~{hpcQueue}"
#         hpcRuntimeMinutes : runtimeMinutes
#         zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
#         preemptible : maxPreemptAttempts
#     }
# }

# task spectra {
#     input {
# Int diskGB = 0
# Int nThreads = 24
# Int gbRAM = 62
# String hpcQueue = "norm"
# Int runtimeMinutes = 240
# String spectraDocker = ""
# Int maxPreemptAttempts = 3
#     }
#     command {
#     }
#     output {
#     }
#     runtime {
#         docker : "~{strspyDocker}"
#         disks : "local-disk ~{auto_diskGB} SSD"
#         cpu : nThreads
#         memory : "~{gbRAM} GB"
#         hpcMemory : gbRAM
#         hpcQueue : "~{hpcQueue}"
#         hpcRuntimeMinutes : runtimeMinutes
#         zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
#         preemptible : maxPreemptAttempts
#     }
# }

task sniffles2 {
    input {
        File inputBAM
        File inputBAI
        File refTarball
        File? snifflesTandemRepeatBed

        Boolean snifflesOutputReadNames = false

        Int diskGB = 0
        Int nThreads = 24
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        String snifflesDocker = "erictdawson/sniffles2"
        Int maxPreemptAttempts = 3
    }
    String localTarball = basename(refTarball)
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB")) + ceil(size(refTarball, "GB") * 3) + 80 else diskGB
    String ref = basename(refTarball, ".tar")
    String outbase = basename(inputBAM, ".bam")

    command {
        mv ~{refTarball} ~{localTarball} && \
        tar xvf ~{localTarball} && \
        sniffles \
        -i ~{inputBAM} \
        --reference ~{ref} \
        ~{if snifflesOutputReadNames then "--output-rnames" else ""} \
        -v ~{outbase}.sniffles2.vcf.gz

    }
    output {
        File outputVCF = "~{outbase}.sniffles2.vcf.gz"
        File outputTBI = "~{outbase}.sniffles2.vcf.gz.tbi"
    }
    runtime {
        docker : "~{snifflesDocker}"
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

workflow AoU_ONT_VariantCalling {
    input {
        File inputBAM
        File inputBAI
        File refTarball

        ## Sniffles Options
        File? snifflesTandemRepeatBed
        Int snifflesThreads = 24
        Int snifflesGBRAM = 64
        Boolean snifflesOutputReadNames = false

        ## Clair3 args
        File? clairTargetsBed

        ## DeepVariant Runtime Args
        String pbPATH = "pbrun"
        File? pbLicenseBin
        String pbDocker = "erictdawson/parabricks-deepvariant"
        Boolean gvcfMode = false
        Int nGPU_DeepVariant = 4
        String gpuModel_DeepVariant = "nvidia-tesla-t4"
        String gpuDriverVersion_DeepVariant = "460.73.01"
        Int nThreads_DeepVariant = 24
        Int gbRAM_DeepVariant = 120
        Int diskGB_DeepVariant = 0
        Int runtimeMinutes_DeepVariant = 600
        String hpcQueue_DeepVariant = "gpu"

        Int maxPreemptAttempts = 3
    }

    call deepvariant{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputRefTarball=refTarball,
            pbLicenseBin=pbLicenseBin,
            pbPATH=pbPATH,
            gvcfMode=gvcfMode,
            nThreads=nThreads_DeepVariant,
            nGPU=nGPU_DeepVariant,
            gpuModel=gpuModel_DeepVariant,
            gpuDriverVersion=gpuDriverVersion_DeepVariant,
            gbRAM=gbRAM_DeepVariant,
            diskGB=diskGB_DeepVariant,
            runtimeMinutes=runtimeMinutes_DeepVariant,
            hpcQueue=hpcQueue_DeepVariant,
            pbDocker=pbDocker,
            maxPreemptAttempts=maxPreemptAttempts 
    }

    call sniffles2 as sniffles{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            refTarball=refTarball,
            snifflesTandemRepeatBed=snifflesTandemRepeatBed,
            nThreads=snifflesThreads,
            snifflesOutputReadNames=snifflesOutputReadNames,
            gbRAM=snifflesGBRAM
    }

    call clair3 as clair{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            refTarball=refTarball,
            targetsBed=clairTargetsBed
    }



    output {
        ## SV calls from Sniffles
        File snifflesVCF = sniffles.outputVCF
        File snifflesTBI = sniffles.outputTBI

        ## SNV calls from Clair3
        File clairPileupVCF = clair.pileupVCF
        File clairFullAlignmentVCF = clair.fullAlignmentVCF
        File clairMergeVCF = clair.mergeVCF

        ## SNV calls from DeepVariant
        File deepvariantVCF = deepvariant.outputVCF
        File deepvariantTBI = deepvariant.outputTBI

        # ## STRspy calls
        # File strspyVCF = strspy.outputVCF
        # File strspyTBI = strspy.outputVCF
    }
}
