# methylDackel

Workflow to run methylDackel, will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq or EM-seq alignments and extract per-base methylation metrics from them. The extract task generates bedGraph files, by default generates only CpG metrics, option can be set to also generate CHH and CHG metrics. Mbias task generates tsv file for methylation bias metrics and a svg graph for visualizing mbias (only for chromosome 1 here).

## Overview

## Dependencies

* [methyldackel 0.6.1](https://github.com/dpryan79/MethylDackel)


## Usage

### Cromwell
```
java -jar cromwell.jar run methylDackel.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`bam`|File|The bam file for methyl analysis
`bai`|File|The index for input bam
`outputFileNamePrefix`|String|Prefix for output files
`reference`|String|The genome reference build


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`doMbias`|Boolean|true|Whether run Mbias or not 


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`methylDackelExtract.doCHH`|Boolean|false|whether enable CHH metrics
`methylDackelExtract.doCHG`|Boolean|false|whether enable CHG metrics
`methylDackelExtract.mergeContext`|Boolean|false|whether merge context in bedgraph
`methylDackelExtract.minimumuQualityPhred`|Int?|None|minimumu sequencing quality phred score
`methylDackelExtract.minimumMAPQ`|Int?|None|minimum MAPQ score
`methylDackelExtract.minDepth`|Int?|None|region with minimum depth needed to be included in analysis
`methylDackelExtract.timeout`|Int|8|The hours until the task is killed
`methylDackelExtract.memory`|Int|16|The GB of memory provided to the task
`methylDackelExtract.threads`|Int|8|The number of threads the task has access to
`extractChromosomes.timeout`|Int|1|The hours until the task is killed
`extractChromosomes.memory`|Int|1|The GB of memory provided to the task
`extractChromosomes.threads`|Int|1|The number of threads the task has access to
`extractChromosomes.modules`|String|"samtools/1.16.1"|The modules that will be loaded
`methylDackelMbias.timeout`|Int|12|The hours until the task is killed
`methylDackelMbias.memory`|Int|8|The GB of memory provided to the task
`methylDackelMbias.threads`|Int|8|The number of threads the task has access to
`concatMbiasFiles.timeout`|Int|1|The hours until the task is killed
`concatMbiasFiles.memory`|Int|1|The GB of memory provided to the task
`concatMbiasFiles.threads`|Int|1|The number of threads the task has access to
`concatMbiasFiles.modules`|String|"pandas/2.1.3"|The modules that will be loaded


### Outputs

Output | Type | Description | Labels
---|---|---|---
`extract_CpGbedgraph`|File|CpGbedGraph output from methylDackelExtract|vidarr_label: extract_CpGbedgraph
`extract_CHGbedgraph`|File?|CHGbedGraph output from methylDackelExtract|vidarr_label: extract_CHGbedgraph
`extract_CHHbedgraph`|File?|CHHbedGraph output from methylDackelExtract|vidarr_label: extract_CHHbedgraph
`mbias_tsv`|File?|mbias tsv output from methylDackelMbias|vidarr_label: mbias_tsv
`mbias_svg`|File?|svg plot files from methylDackelMbias|vidarr_label: mbias_svg


## Commands
 This section lists command(s) run by methylDackel workflow
 
 * Running methylDackel
 
 ```
         samtools view -H ~{bam} | grep @SQ | cut -f2 | sed 's/SN://' | grep -E -v '(_random|chrUn|chrM|MT|_alt|_fix|_decoy|_PATCH|_HSCHR|NC_|_EBV|EBV|phiX|pUC19|lambda|_scaffold)'
 ```
 ```
         set -euo pipefail
         MethylDackel extract ~{filterMAPQ} ~{filterQalityPhred} ~{filterminDepth} ~{optionMergeContext} ~{optionCHH} ~{optionCHG} -@ ~{threads} ~{fasta} ~{bam} -o ~{outputFileNamePrefix}.methyldackel
         
 ```
 ```
         MethylDackel mbias --txt -r ~{chr} ~{fasta} ~{bam} ~{outputFileNamePrefix}.mbias > output_mbias.tsv
         tar  -czf ~{outputFileNamePrefix}_mbias.svgs.tar.gz *.svg
 ```
 ```
         python3<<CODE
 
         import sys
         import pandas as pd
 
         dfs = []
         input_files = ['~{sep="', '" select_all(inputTsvs)}']
         columns = ['Strand', 'Read', 'Position', 'nMethylated', 'nUnmethylated']
         for file in input_files:
             df = pd.read_csv(file, sep='\t', skiprows=1, names=columns)  # Skip header
             dfs.append(df)
 
         combined_df = pd.concat(dfs, ignore_index=True)
 
         # Group by Strand, Read, and Position, and sum the methylation counts
         aggregated_df = combined_df.groupby(['Strand', 'Read', 'Position'], as_index=False).agg({
             'nMethylated': 'sum',
             'nUnmethylated': 'sum'
         }).sort_values(['Strand', 'Read', 'Position'])
 
         with open("~{outputFileNamePrefix}.mbias.tsv", 'w') as f:
             aggregated_df.to_csv(f, sep='\t', index=False)
         CODE
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
