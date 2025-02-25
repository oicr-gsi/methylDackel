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

    call methylDackel_extract {
        input:
            bam = bam,
            bai = bai,
            prefix = outputFileNamePrefix,
            fasta = ref.fasta,
            modules = "methyldackel/0.6.1 ~{ref.genomeModule}"
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
            description: "bedGraph output from methyDackel extract",
            vidarr_label: "extract_bedgraph"
        }
        }
    }

    output {
        File extract_bedgraph = methylDackel_extract.out
    }
}

task methylDackel_extract {
    input {
        File bam
        File bai
        String prefix
        String fasta
        Int timeout = 6
        Int memory = 8
        Int threads = 8
        String modules
    }

    parameter_meta {
        bam: "The bam file to analyze"
        bai: "The .bai index of the bam file"
        prefix: "File prefix"
        fasta: "FastA file used for alignment"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        MethylDackel extract --mergeContext -@ ~{threads} ~{fasta} ~{bam} -o ~{prefix}.methyldackel
        gzip ~{prefix}.methyldackel_CpG.bedGraph
    >>>

    output {
        File out = "~{prefix}.methyldackel_CpG.bedGraph.gz"
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

