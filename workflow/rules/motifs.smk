motifs = {"CTCF": "MA0139.1", "ATF3": "MA0605.2", "MAFK": "MA0496.1"}
jaspar_address = "http://jaspar.genereg.net/api/v1/matrix/"

rule get_meme:
    output:
        "results/motifs/{condition}.meme"
    run:
        shell('wget {jaspar_address}%s.meme -O {output} --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"' % motifs[wildcards.condition])


rule sort_peaks:
    input:
        "results/{caller}/{regions}/{filename}"
    output:
        "results/{caller}/sorted_{regions}/{filename}"
    shell:
        "sort -nr -k5 {input} > {output}"

rule get_peak_sequences:
    input:
        peaks="results/{caller}/{regions}/{sample}.narrowPeak",
    output:
        "results/{caller}/{regions}_fasta/{sample}.fa"
    params:
        reference=config["reference_fasta"]
    shell:
        "bedtools getfasta -fi {params.reference} -bed {input.peaks} > {output}"

rule motif_enrichment:
    input:
        fasta="results/{caller}/{regions}_fasta/{celltype}_{condition}.fa",
        meme="results/motifs/{condition}.meme"
    output:
        multiext("results/{caller}/motif_matches_{regions}/{celltype}_{condition}/fimo", ".html", ".xml", ".tsv", ".gff") 
    shell:
        "fimo --oc results/{wildcards.caller}/motif_matches_{wildcards.regions}/{wildcards.celltype}_{wildcards.condition}/ {input.meme} {input.fasta}"

rule motif_plot:
    input:
        matches="results/{caller}/motif_matches_{regions}/{sample}/fimo.tsv",
        peaks="results/{caller}/sorted_{regions}/{sample}.narrowPeak"
    output:
        "results/{caller}/motif_plots/{regions}/{sample}.png",
        "results/{caller}/motif_plots/{regions}/{sample}.pkl"
    script:
        "../scripts/motifplot.py"

rule join_plot:
    input:
        expand("results/{caller}/motif_plots/{{regions}}/{{sample}}.pkl", caller=config["callers"])
    output:
        report("results/plots/motif/{regions}/{sample}.png", category="Motif plots")
    script:
        "../scripts/joinplots.py"

