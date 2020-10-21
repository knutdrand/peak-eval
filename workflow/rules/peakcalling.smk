import pandas as pd
import os
macs_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.narrowPeak"]


rule macs_wo_control:
    input:
        treatment=lambda w: [f"results/reads/{sample}.bed.gz" for sample in samples.index[(samples["condition"]==w.condition) & (samples["celltype"]==w.celltype)]]
    output:
        expand("results/macs2_folder/{{celltype}}_{{condition}}_{filetype}", filetype=macs_output)
    shell:
        """macs2 callpeak -t {input.treatment} -g 3000000 --bdg --outdir results/macs2_folder/ -n {wildcards.celltype}_{wildcards.condition} --nomodel --extsize 176"""

rule hmmcaller_wo_control:
    input:
        treatment="results/macs2_folder/{sample}_treat_pileup.bdg",
        control="results/macs2_folder/{sample}_control_lambda.bdg"
    output:
        "results/hmmcaller/peaks/{sample}.narrowPeak"
    shell:
        "hmmcaller {input.treatment} {input.control} {output}"

rule filter_hmmcaller:
    input:
        "results/hmmcaller/peaks/{sample}.narrowPeak"
    output:
        "results/filtered_hmmcaller/peaks/{sample}.narrowPeak"
    shell:
        "awk '{{if ($3-$2>=176) print}}' {input} > {output}"

rule callpeak_w_control:
    input:
        treatment=lambda w: expand("results/reads/{sample}.bed.gz", sample==samples.index[samples["condition"]==w.condition & samples["celltype"]==w.celltype]),
        control=lambda w: expand("results/reads/{sample}.bed.gz", 
                                 sample==samples.index[samples["condition"]=="input" & samples["celltype"]==w.celltype])
    output:
        expand("results/controlled_macs2_folder/{{celltype}}_{{condition}}_{filetype}", filetype=macs_output)
    shell:
        """macs2 callpeak -t {input.treatment} -c {input.control} -g hg --bdg --outdir macs2_folder/ -n {wildcards.celltype}_{wildcards.condition}"""

rule move_macs_peaks:
    input: 
        "results/macs2_folder/{sample}_peaks.narrowPeak"
    output:
        "results/macs2/peaks/{sample}.narrowPeak"
    shell:
        "mv {input} {output}"

rule encode_download:
    output:
        "results/reads/{sample}.bam"
    shell:
        temp("wget https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bam -O {output}")

rule region_bamtobed:
    input:
        "results/reads/{sample}.bam",
    output:
        temp("results/unsorted_reads/{sample}.bed")
    shell:
        """bedtools bamtobed -i {input} | awk '{{if ($1=="chr1" && $3<3000000) print}}' > {output}"""


rule sortbed:
    input:
        "results/unsorted_reads/{sample}.bed"
    output:
        "results/reads/{sample}.bed.gz"
    shell:
        "sort -k 1,1 -k2,2n {input} | gzip > {output}"
