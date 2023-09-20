version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/vcf.wdl" as vcf
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/tasks/deepvariant.wdl" as deepvariant

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
        
        String clairDocker = "erictdawson/clair3:latest"
        RuntimeAttributes runtime_attributes
        GPUAttributes gpu_attributes

    }
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + ceil(size(refTarball, "GB")) + ceil(size(inputBAI, "GB")) + 65 else runtime_attributes.diskGB
    String ref = basename(refTarball, ".tar")
    String outbase = basename(inputBAM, ".bam")
    command {
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        tar xvf ~{refTarball}

        /opt/bin/run_clair3.sh ~{"--sample_name " + sampleName} --ref_fn=~{ref} --threads=~{runtime_attributes.nThreads} --platform=~{platform} --model_path=~{modelPath} --output=~{outbase}.clair3 --bam_fn=~{inputBAM} ~{if gvcfMode then "--gvcf" else ""} ~{if phaseMode then "--enable_phasing" else ""} ~{"--bed_fn=" + targetsBed}

    }
    output {
        File pileupVCF = "~{outbase}.clair3/pileup.vcf.gz"
        File fullAlignmentVCF = "~{outbase}.clair3/full_alignment.vcf.gz"
        File mergeVCF = "~{outbase}.clair3/merge_output.vcf.gz"
    }
    runtime {
        docker : "~{clairDocker}"
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

task mosdepth {
    input {
        File inputBAM
        File inputBAI
        Int windowSize = 20
        Int minMAPQ = 20
        Boolean noPerBase = true
        Boolean useMedian = false
        String? thresholds
        String? quantize
        
        String mosdepthDocker = "erictdawson/mosdepth"
        RuntimeAttributes runtime_attributes = {
            "diskGB": 0,
            "nThreads": 4,
            "gbRAM": 11,
            "hpcQueue": "norm",
            "runtimeMinutes": 600,
            "maxPreemptAttempts": 3
        }
    }
    String outbase = basename(inputBAM, ".bam")
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputBAM, "GB") * 3.2) + 80 else runtime_attributes.diskGB

    command {
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mosdepth \
            ~{"--by " + windowSize} \
            --threads ~{runtime_attributes.nThreads} \
            ~{if noPerBase then "--no-per-base" else ""} \
            ~{"--mapq " + minMAPQ} \
            ~{"--thresholds " + thresholds} \
            ~{"--quantize " + quantize} \
            ~{if useMedian then "--use-median" else ""} \
            ~{outbase} \
            ~{inputBAM} && \
            mkdir mosdepth && \
            mv ~{outbase}.mosdepth* mosdepth/ && \
            tar cvf mosdepth.tar mosdepth

    }
    output {
        File outputTarball = "mosdepth.tar"
        # File outputRegionDepthFile = "~{outbase}.mosdepth.region.dist.txt"
        # File outputGlobalDepthFile = "~{outbase}.mosdepth.global.dist.txt"
        # File outputSummaryFile = "~{outbase}.mosdepth.summary.txt"
    }
    runtime {
        docker : "~{mosdepthDocker}"
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


task spectre {
    input {
        File mosDepthTarball
        File inputRefTarball
        File sampleName
        Int binSize

        File? snvVCF
        File? blacklistBED
        Int? ploidy
        Int? nSize
        Boolean populationMode = false
        Int? minCNVLen

        String spectreDocker = "erictdawson/spectre"
        RuntimeAttributes runtime_attributes = {
            "diskGB": 0,
            "nThreads": 4,
            "gbRAM": 11,
            "hpcQueue": "norm",
            "runtimeMinutes": 600,
            "maxPreemptAttempts": 3
        }
    }

    String ref = basename(inputRefTarball, ".tar")
    Int auto_diskGB = 0
    String outbase = sampleName + ".spectre"
    String mosdepthCoverageDirectory = ""
    command {
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        mv ~{mosDepthTarball} . && \
        tar xvf ~{mosDepthTarball} && \
        mv ~{inputRefTarball} . && \
        tar xvf ~{inputRefTarball} && \
        python spectre/spectre.py CNVCaller \
        --bin-size ~{binSize} \
        --reference ~{ref} \
        --sample-id sampleName \
        --coverage ~{mosdepthCoverageDirectory} \
        --output-dir ~{outbase} \
        ~{"--threads " + runtime_attributes.nThreads} \
        ~{"--blacklist " + blacklistBED} \
        ~{"--ploidy " + ploidy} \
        ~{"--snv " + snvVCF} \
        ~{"--n-size " + nSize} \
        ~{"--min_cnv_len " + minCNVLen} \
        ~{if populationMode then "--population" else ""} && \
        tar cvf ~{outbase}.tar ~{outbase}
    }
    output {
        File outputTarball = "~{outbase}.tar"
    }
    runtime {
        docker : "~{spectreDocker}"
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

task sniffles2 {
    input {
        File inputBAM
        File inputBAI
        File refTarball
        String sampleName
        File? tandemRepeatBed

        Boolean mosaicMode = false
        Boolean outputReadNames = false
        Int minSVLen = 50

        String snifflesDocker = "erictdawson/sniffles2"
        RuntimeAttributes runtime_attributes = {
            "diskGB": 0,
            "nThreads": 24,
            "gbRAM": 62,
            "hpcQueue": "norm",
            "runtimeMinutes": 240,
            "maxPreemptAttempts": 3
        }
    }
    String localTarball = basename(refTarball)
    Int auto_diskGB = if runtime_attributes.diskGB == 0 then ceil(size(inputBAM, "GB")) + ceil(size(refTarball, "GB") * 3) + 80 else runtime_attributes.diskGB
    String ref = basename(refTarball, ".tar")
    String outbase = basename(inputBAM, ".bam")

    command {
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        mv ~{refTarball} ~{localTarball} && \
        tar xvf ~{localTarball} && \
        sniffles \
        ~{"--minsvlen " + minSVLen} \
        --threads ~{runtime_attributes.nThreads} \
        --reference ~{ref} \
        --input ~{inputBAM} \
        --sample-id ~{sampleName} \
        --vcf ~{outbase}.sniffles2.vcf.gz \
        --snf ~{outbase}.sniffles2.snf
        ~{if mosaicMode then "--mosaic" else ""} \
        ~{"--tandem-repeats " + tandemRepeatBed} \
        ~{if outputReadNames then "--output-rnames" else ""} 

    }
    output {
        File outputVCF = "~{outbase}.sniffles2.vcf.gz"
        File outputTBI = "~{outbase}.sniffles2.vcf.gz.tbi"
        File outputSNF = "~{outbase}.sniffles2.snf"
    }
    runtime {
        docker : "~{snifflesDocker}"
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
        String pbDocker = "nvcr.io/nv-parabricks-dev/clara-parabricks:4.2.0-1.beta3"
        Boolean gvcfMode = false
        File? deepvariantModelFile
        String deepvariantMode = "ont"

        ## MosDepth inputs
        Int windowSize = 20
        Int minMAPQ = 20

        Int maxPreemptAttempts = 3
    }

    RuntimeAttributes deepvariant_attributes = {
        "diskGB": 0,
        "nThreads": 48,
        "gbRAM": 160,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    GPUAttributes deepvariant_gpu_attributes = {
        "gpuModel" : "nvidia-tesla-t4",
        "nGPU" : 4,
        "gpuDriverVersion": "535.104.05"
    }

    call deepvariant.deepvariant as deepvariant{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            inputRefTarball=refTarball,
            gvcfMode=gvcfMode,
            runtime_attributes=deepvariant_attributes,
            gpu_attributes=deepvariant_gpu_attributes
    }

    call vcf.compressAndIndexVCF as compress_deepVariant {
        input:
            inputVCF=deepvariant.outputVCF
    }

    RuntimeAttributes sniffles_attributes = {
        "diskGB": 0,
        "nThreads": 24,
        "gbRAM": 62,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    call sniffles2 as sniffles{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            sampleName=sampleName,
            refTarball=refTarball,
            tandemRepeatBed=snifflesTandemRepeatBed,
            outputReadNames=snifflesOutputReadNames,
            runtime_attributes=sniffles_attributes
    }

    RuntimeAttributes clair_attributes = {
        "diskGB": 0,
        "nThreads": 24,
        "gbRAM": 87,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "maxPreemptAttempts": 3,
    }

    GPUAttributes clair_gpu_attributes = {
        "gpuModel" : "nvidia-tesla-t4",
        "nGPU" : 4,
        "gpuDriverVersion": "535.104.05"
    }

    call clair3 as clair{
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            refTarball=refTarball,
            targetsBed=clairTargetsBed,
            runtime_attributes=clair_attributes
    }


    call mosdepth {
        input:
            inputBAM=inputBAM,
            inputBAI=inputBAI,
            windowSize=windowSize,
            minMAPQ=minMAPQ
    }

    call spectre {
        input:
            mosDepthTarball=mosdepth.outputTarball,
            inputRefTarball=refTarball,
            sampleName=sampleName,
            binSize=windowSize
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
        File deepvariantVCF = compress_deepVariant.outputVCFGZ
        File deepvariantTBI = compress_deepVariant.outputTBI

        # Spectre CNV calls
        File spectreOutput = spectre.outputTarball

    }
}
