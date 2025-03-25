version 1.0

struct GenomeResources {
    String fasta
    String genomeModule
}

workflow methylDackel {
    input {
        File bam
        File bai
        String outputFileNamePrefix
        String reference
        Boolean doMbias = true
    }

    parameter_meta {
        bam: "The bam file for methyl analysis"
        bai: "The index for input bam"
        outputFileNamePrefix: "Prefix for output files"
        reference: "The genome reference build"
    }

    Map[String, GenomeResources] resources = {
        "hg38": {
            "fasta": "$HG38_EM_SEQ_ROOT/hg38_random.fa",
            "genomeModule": "hg38-em-seq/p12-2022-10-17"
        }
    }

    GenomeResources ref = resources[reference]

    call methylDackelExtract {
        input:
            bam = bam,
            bai = bai,
            outputFileNamePrefix = outputFileNamePrefix,
            fasta = ref.fasta,
            modules = "methyldackel/0.6.1 ~{ref.genomeModule}"
    }

    if ( doMbias ){
        call extractChromosomes{
            input:
            bam = bam
        }

        scatter ( chr in extractChromosomes.chromosomes ) {
            call methylDackelMbias {
                input:
                    bam = bam,
                    bai = bai,
                    chr = chr,
                    fasta = ref.fasta,
                    modules = "methyldackel/0.6.1 samtools/1.16.1 ~{ref.genomeModule}"
            }
        }      
    
        Array[File?] mbiasTsvs = select_first([methylDackelMbias.mbias_tsv])
        Array[Array[File?]] mbiasSvg = select_first([methylDackelMbias.mbias_svg_files])

        call concatenateTsvFiles {
            input:
                inputTsvs = mbiasTsvs,
                outputFileNamePrefix = outputFileNamePrefix
        }
    }

    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "Workflow to run bwa-meth, the fast aligner for EM-seq/BS-Seq reads. Prior to alignment, adatper trimming and quality filtering are performed. Readgroup information to be injected into the bam header needs to be provided.  The workflow can also split the input data into a requested number of chunks, align each separately then merge the separate alignments into a single bam file.  This decreases the workflow run time. Final bam file also applied markDuplicates."
        dependencies: [
         {
            name: "methyldackel/0.6.1",
            url: "https://github.com/dpryan79/MethylDackel"
                        }
        ]
        output_meta: {
            extract_bedgraph: {
                description: "bedGraph output from methylDackelExtract",
                vidarr_label: "extract_bedgraph"
            },
            combined_mbias_tsv: {
                description: "mbias tsv output from methylDackelMbias",
                vidarr_label: "combined_mbias_tsv"
            },
            mbias_svg: {
                description: "svg plot files from methylDackelMbias",
                vidarr_label: "mbias_svg_files"
            }
        }
    }

    output {
        Array[File] extract_bedgraph = methylDackelExtract.out
        File? mbias_tsv = concatenateTsvFiles.combinedTsv
        Array[File?] mbias_svg = mbiasSvg[0]
    }
}

task extractChromosomes {
    input {
        File bam
        Int timeout = 1
        Int memory = 1
        Int threads = 1
        String modules = "samtools/1.16.1"
    }

    parameter_meta {
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        samtools view -H ~{bam} | grep @SQ | cut -f2 | sed 's/SN://' | grep -E -v '(_random|chrUn|chrM|MT|_alt|_fix|_decoy|_PATCH|_HSCHR|NC_|_EBV|phiX|pUC19|lambda|_scaffold)'
    >>>

    output {
        Array[String] chromosomes = read_lines(stdout())
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task methylDackelExtract {
    input {
        File bam
        File bai
        String outputFileNamePrefix
        String fasta
        Boolean doCHH = false
        Boolean doCHG = false
        Boolean mergeContext = false
        Int minimumuQalityPhred = 5
        Int minimumMAPQ = 10
        Int timeout = 12
        Int memory = 8
        Int threads = 8
        String modules
    }

    parameter_meta {
        bam: "The bam file to analyze"
        bai: "The .bai index of the bam file"
        outputFileNamePrefix: "Output file name prefix"
        fasta: "FastA file used for alignment"
        doCHH: "whether enable CHH metrics"
        doCHG: "whether enable CHG metrics"
        mergeContext: "whther merge context in bedgraph"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }
    String optionCHH = if doCHH then "--CHH" else ""
    String optionCHG = if doCHG then "--CHG" else ""
    String optionMergeContext = if mergeContext then "--mergeContext" else ""

    command <<<
        set -euo pipefail
        MethylDackel extract -q {minimumMAPQ} -p ~{minimumuQalityPhred} ~{optionMergeContext} ~{optionCHH} ~{optionCHG} -@ ~{threads} ~{fasta} ~{bam} -o ~{outputFileNamePrefix}.methyldackel
        gzip *.bedGraph
        mkdir -p ~{outputFileNamePrefix}_extract_bedGraph
        mv *.bedGraph.gz ~{outputFileNamePrefix}_extract_bedGraph
    >>>

    output {
        Array[File] out = glob("~{outputFileNamePrefix}_extract_bedGraph/*.bedGraph.gz")
    }

    meta {
        output_meta: {
            out: "The compressed MethylDackel result bedGraph"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task methylDackelMbias {
    input {
        File bam
        File bai
        String chr
        String fasta
        String modules
        Int timeout = 48
        Int memory = 8
        Int threads = 8
    }

    parameter_meta {
        bam: "The bam file to analyze"
        bai: "The .bai index of the bam file"
        chr: "The region to call methylDackel mbias"
        fasta: "FastA file used for alignment"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
       MethylDackel mbias --txt -r ~{chr} ~{fasta} ~{bam} output.mbias > output_mbias.tsv
    >>>

    output {
        File? mbias_tsv = "output_mbias.tsv"
        Array[File?] mbias_svg_files = ["output.mbias_OT.svg", "output.mbias_OB.svg"]
    }
    
    meta {
        output_meta: {
            mbias_tsv: "mbias tsv output from methylDackelMbias",
            mbias_svg_files: "svg plot files from methylDackelMbias"
        }
    }
    
    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task concatenateTsvFiles {
    input {
        Array[File?] inputTsvs
        String outputFileNamePrefix
    }

    command <<<
        for tsv in ~{sep=' ' select_all(inputTsvs)}; do
            cat "$tsv"
        done >> ~{outputFileNamePrefix}.combined.tsv
    >>>

    output {
        File combinedTsv = "~{outputFileNamePrefix}.combined.tsv"
    }
}