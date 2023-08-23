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
        String clairDocker = "erictdawson/clair3:latest"
        Int diskGB = 0
        Int nGPU = 2
        String gpuModel = "nvidia-tesla-t4"
        Int maxPreemptAttempts = 3
        Int runtimeMinutes = 300
        String hpcQueue = "norm"
        Int gbRAM = 87

    }
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(refTarball, "GB")) + ceil(size(inputBAI, "GB")) + 65 else diskGB
    String ref = basename(refTarball, ".tar")
    String outbase = basename(inputBAM, ".bam")
    command {
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        tar xvf ~{refTarball}

        # REF=$(basename ~{ref})
        # REF_IDX=$(basename ~{ref}.fai)
        # ln -s ~{ref} ./$REF
        # ln -s ~{ref}.fai ./$REF_IDX

        # BAM=$(basename ~{inputBAM})
        # BAI=$(basename ~{inputBAI})
        # ln -s ~{inputBAM} ./$BAM
        # ln -s ~{inputBAI} ./$BAI

        /opt/bin/run_clair3.sh --ref_fn=~{ref} --threads=~{nThreads} --platform=~{platform} --model_path=~{modelPath} --output=~{outbase}.clair3 --bam_fn=~{inputBAM} 

    # ~{"--bed_fn=" + targetsBed} \
    # ~{"--sample_name=" + sampleName} \
    # ~{if gvcfMode then "--gvcf" else ""} \
    # ~{if phaseMode then "--enable_phasing" else ""} \
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
        gpuType : "~{gpuModel}"
        gpuCount : nGPU
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
        File? modelFile
        String pbPATH = "pbrun"
        File? pbLicenseBin
        String pbDocker = "nvcr.io/nvidia/clara/clara-parabricks:4.1.0-1"
        Boolean gvcfMode = false
        Int nGPU = 4
        String gpuModel = "nvidia-tesla-t4"
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
        nvidia-smi && \
        time /usr/local/parabricks/binaries/bin/deepvariant \
        ~{ref} \
        ~{inputBAM} 2 2 \
        -o ~{outVCF} \
        -n 6 \
        --model ~{modelFile} \
        -long_reads \
        --sort_by_haplotypes \
        --add_hp_channel \
        --parse_sam_aux_fields \
        -norealign_reads \
        --vsc_min_fraction_indels 0.06 \
        --track_ref_reads \
        --phase_reads \
        --pileup_image_width 199 \
        --max_reads_per_partition 600 \
        --partition_size 25000 \
        --vsc_min_count_snps 2 \
        --vsc_min_count_indels 2 \
        --vsc_min_fraction_snps 0.12 \
        --min_mapping_quality 1 \
        --min_base_quality 10 \
        --alt_aligned_pileup diff_channels \
        --variant_caller VERY_SENSITIVE_CALLER \
        --dbg_min_base_quality 15 \
        --ws_min_windows_distance 80 \
        --aux_fields_to_keep HP \
        --p_error 0.001 \
        --max_ins_size 10 && \
        bgzip ~{outVCF} && \
        tabix ~{outVCF}.gz




        # time ${pbPATH} deepvariant \
        # --x3 \
        # ~{if gvcfMode then "--gvcf " else ""} \
        # --ref ${ref} \
        # --in-bam ${inputBAM} \
        # --out-variants ~{outVCF} \
        # ~{"--pb-model-file " + modelFile} \
        # --mode pacbio \
        # ~{"--license-file " + pbLicenseBin} && \
        # bgzip ~{outVCF} && \
        # tabix ~{outVCF}.gz
    }
    output {
        File outputVCF = "~{outVCF}.gz"
        File outputTBI = "~{outVCF}.gz.tbi"
    }
    runtime {
        docker : "~{pbDocker}"
        disks : "local-disk ~{auto_diskGB} SSD"
        cpu : nThreads
        memory : "~{gbRAM} GB"
        hpcMemory : gbRAM
        hpcQueue : "~{hpcQueue}"
        hpcRuntimeMinutes : runtimeMinutes
        gpuType : "~{gpuModel}"
        gpuCount : nGPU
        zones : ["us-central1-a", "us-central1-b", "us-central1-c"]
        preemptible : maxPreemptAttempts
    }
}

task mosdepth {
    input {
        File inputBAM
        File inputBAI
        Int windowSize = 20
        Boolean noPerBase = true

        Int diskGB = 0
        Int nThreads = 12
        Int gbRAM = 62
        String hpcQueue = "norm"
        Int runtimeMinutes = 240
        String mosdepthDocker = "erictdawson/mosdepth"
        Int maxPreemptAttempts = 3
    }
    String outbase = basename(inputBAM, ".bam")
    Int auto_diskGB = if diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + 80 else diskGB

    command {
        mosdepth \
            ~{"--by " + windowSize} \
            --threads ~{nThreads} \
            ~{if noPerBase then "--no-per-base" else ""} \
            --mapq 20 \
            ~{outbase} \
            ~{inputBAM}

    }
    output {
        File outputRegionDepthFile = "~{outbase}.mosdepth.region.dist.txt"
        File outputGlobalDepthFile = "~{outbase}.mosdepth.global.dist.txt"
        # File outputSummaryFile = "~{outbase}.mosdepth.summary.txt"
    }
    runtime {
        docker : "~{mosdepthDocker}"
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


# task spectre {
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
        String sampleName
        File? snifflesTandemRepeatBed

        Boolean mosaicMode = false
        Boolean snifflesOutputReadNames = false
        Int minSVLen = 50

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
        --minsvlen ~{minSVLen} \
        --threads ~{nThreads} \
        --reference ~{ref} \
        --input ~{inputBAM} \
        ~{if mosaicMode then "--mosaic" else ""} \
        --sample-id ~{sampleName} \
        ~{if snifflesOutputReadNames then "--output-rnames" else ""} \
        --vcf ~{outbase}.sniffles2.vcf.gz \
        --snf ~{outbase}.sniffles2.snf

    }
    output {
        File outputVCF = "~{outbase}.sniffles2.vcf.gz"
        File outputTBI = "~{outbase}.sniffles2.vcf.gz.tbi"
        File outputSNF = "~{outbase}.sniffles2.snf"
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
	    String sampleName
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
        String pbDocker = "nvcr.io/nvidia/clara/clara-parabricks:4.1.0-1"
        Boolean gvcfMode = false
        Int nGPU_DeepVariant = 4
        String gpuModel_DeepVariant = "nvidia-tesla-t4"
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
            sampleName=sampleName,
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


    call mosdepth {
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI
    }




    output {
        ## SV calls from Sniffles
        File snifflesVCF = sniffles.outputVCF
        File snifflesTBI = sniffles.outputTBI

        ## SNV calls from Clair3
        File clairPileupVCF = clair.pileupVCF
        File clairFullAlignmentVCF = clair.fullAlignmentVCF
        File clairMergeVCF = clair.mergeVCF

        # SNV calls from DeepVariant
        File deepvariantVCF = deepvariant.outputVCF
        File deepvariantTBI = deepvariant.outputTBI

        # ## STRspy calls
        # File strspyVCF = strspy.outputVCF
        # File strspyTBI = strspy.outputVCF
    }
}
