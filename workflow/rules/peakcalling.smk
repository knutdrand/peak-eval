import pandas as pd
import os
macs_output=["treat_pileup.bdg", "control_lambda.bdg", "peaks.narrowPeak"]
genomesize = config.get("max_pos", "hs")

rule macs_wo_control:
    input:
        treatment=lambda w: [f"results/reads/{sample}.bed.gz" for sample in samples.index[(samples["condition"]==w.condition) & (samples["celltype"]==w.celltype)]]
    output:
        expand("results/macs2_folder/{{celltype}}_{{condition}}_{filetype}", filetype=macs_output)
    shell:
        """macs2 callpeak -t {input.treatment} -g {genomesize} --bdg --outdir results/macs2_folder/ -n {wildcards.celltype}_{wildcards.condition} --nomodel --extsize 176"""

rule hmmacs_wo_control:
    input:
        treatment="results/macs2_folder/{sample}_treat_pileup.bdg",
        control="results/macs2_folder/{sample}_control_lambda.bdg"
    output:
        "results/hmmacs/peaks/{sample}.narrowPeak"
    shell:
        "hmmacs {input.treatment} {input.control} {output}"

rule filter_hmmacs:
    input:
        "results/hmmacs/peaks/{sample}.narrowPeak"
    output:
        "results/filtered_hmmacs/peaks/{sample}.narrowPeak"
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
        """macs2 callpeak -t {input.treatment} -c {input.control} -g {genomesize} --bdg --outdir macs2_folder/ -n {wildcards.celltype}_{wildcards.condition}"""

rule move_macs_peaks:
    input: 
        "results/macs2_folder/{sample}_peaks.narrowPeak"
    output:
        "results/macs2/peaks/{sample}.narrowPeak"
    shell:
        "mv {input} {output}"

rule summit_regions:
    input:
        "results/{caller}/peaks/{sample}.narrowPeak"
    output:
        "results/{caller}/summit_regions/{sample}.narrowPeak"
    shell:
        """awk '{{OFS="\t"}}{{print $1, $2+$10-200, $2+$10+200, ".", $5}}' {input} > {output}"""
