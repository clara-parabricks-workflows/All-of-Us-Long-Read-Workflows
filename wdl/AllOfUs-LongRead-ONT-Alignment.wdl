version 1.0

import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/util/attributes.wdl"
import "https://raw.githubusercontent.com/clara-parabricks-workflows/parabricks-wdl/long-read/wdl/long-read/minimap2.wdl" as alignFASTQ

workflow AoU_ONT_Alignment {
    input {
        File inputFASTQ
        File inputReference
        File? referenceIndex
        String sampleName
        String mm2Preset = "map-ont"
        String mm2Flags = " -Y "
        Boolean addMDTag = true
        Int mapThreads = 32
        Int sortThreads = 6
        Int minimap_RAM = 120
        Int sortRAM_per_thread = 6
        String minimapDocker = "erictdawson/minimap2"
        String minimap_queue = "norm"
        Int minimap_runtime_max = 600
    }

        RuntimeAttributes attributes = {
        "diskGB": 0,
        "nThreads": mapThreads + sortThreads + 2,
        "gbRAM": minimap_RAM + sortThreads*sortRAM_per_thread + 24,
        "hpcQueue": "norm",
        "runtimeMinutes": 600,
        "gpuDriverVersion": "535.104.05",
        "maxPreemptAttempts": 3,
    }

    call alignFASTQ.minimap2 as minimap2{
        input:
            inputFASTQ=inputFASTQ,
            inputReference=inputReference,
            referenceIndex=referenceIndex,
            sampleName=sampleName,
            addMDTag=addMDTag,
            mm2Flags=mm2Flags,
            mm2Preset=mm2Preset,
            mapThreads=mapThreads,
            sortThreads=sortThreads,
            sortRAM_per_thread=sortRAM_per_thread,
            minimapDocker=minimapDocker,
            runtime_attributes=attributes
    }


    output {
        File alignedBAM = minimap2.outputBAM
        File alignedBAI = minimap2.outputBAI
    }
}
