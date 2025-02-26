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

    call methylDackelExtract {
        input:
            bam = bam,
            bai = bai,
            outputFileNamePrefix = outputFileNamePrefix,
            fasta = ref.fasta,
            modules = "methyldackel/0.6.1 ~{ref.genomeModule}"
    }

    call methylDackelMbias {
        input:
            bam = bam,
            bai = bai,
            outputFileNamePrefix = outputFileNamePrefix,
            fasta = ref.fasta,
            modules = "methyldackel/0.6.1 samtools/1.16.1 ~{ref.genomeModule}"
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
            mbias_svg_files: {
                description: "svg plot files from methylDackelMbias",
                vidarr_label: "mbias_svg_files"
            }
        }
    }

    output {
        Array[File] extract_bedgraph = methylDackelExtract.out
        File combined_mbias_tsv = methylDackelMbias.combined_mbias_tsv
        Array[File] mbias_svg_files = methylDackelMbias.mbias_svg_files
    }
}

task methylDackelExtract {
    input {
        File bam
        File bai
        String outputFileNamePrefix
        String fasta
        Int timeout = 6
        Int memory = 8
        Int threads = 8
        String modules
    }

    parameter_meta {
        bam: "The bam file to analyze"
        bai: "The .bai index of the bam file"
        outputFileNamePrefix: "Output file name prefix"
        fasta: "FastA file used for alignment"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        MethylDackel extract --mergeContext --CHH --CHG -@ ~{threads} ~{fasta} ~{bam} -o ~{outputFileNamePrefix}.methyldackel
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
        String outputFileNamePrefix
        File bam
        File bai
        String fasta
        String modules
        Int timeout = 6
        Int memory = 8
        Int threads = 8
    }

    parameter_meta {
        bam: "The bam file to analyze"
        bai: "The .bai index of the bam file"
        outputFileNamePrefix: "Output file name prefix"
        fasta: "FastA file used for alignment"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > ~{outputFileNamePrefix}_combined_mbias.tsv

        chrs=(`samtools view -H ~{bam}| grep @SQ | cut -f 2 | sed 's/SN://' | grep -v _random | grep -v chrUn | sed 's/|/\\|/'`)

        for chr in ${chrs[*]}; do
            for context in CHH CHG CpG; do
                arg=''
                if [ $context = 'CHH' ]; then
                    arg='--CHH --noCpG'
                elif [ $context = 'CHG' ]; then
                    arg='--CHG --noCpG'
                fi

                join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
                <( \
                    MethylDackel mbias --noSVG $arg -r $chr ~{fasta} ~{bam}| \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                <( \
                    MethylDackel mbias --noSVG --keepDupes -F 2816 $arg -r $chr ~{fasta} ~{bam}| \
                    tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
                ) \
                | sed "s/^/${chr}\t${context}\t/" \
                >> ~{outputFileNamePrefix}_combined_mbias.tsv
            done
        done

        # Generate SVG files for trimming checks
        MethylDackel mbias --noCpG --CHH --CHG -r ${chrs[0]} ~{fasta} ~{bam} ~{outputFileNamePrefix}.mbias_chn
        for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" $f; done;

        MethylDackel mbias -r ${chrs[0]} ~{fasta} ~{bam} ~{outputFileNamePrefix}.mbias_cpg
        for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" $f; done;
        mkdir -p ~{outputFileNamePrefix}_mbias_svg_files
        mv *.svg ~{outputFileNamePrefix}_mbias_svg_files
    >>>

    output {
        File combined_mbias_tsv = "~{outputFileNamePrefix}_combined_mbias.tsv"
        Array[File] mbias_svg_files = glob("~{outputFileNamePrefix}_mbias_svg_files/*.svg")
    }
    
    meta {
        output_meta: {
            combined_mbias_tsv: "mbias tsv output from methylDackelMbias",
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
