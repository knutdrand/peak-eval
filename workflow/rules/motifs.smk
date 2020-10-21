motifs = {"CTCF": "MA0139.1"}
jaspar_address = "http://jaspar.genereg.net/api/v1/matrix/"

rule get_meme:
    output:
        "results/motifs/{condition}.meme"
    run:
        shell('wget {jaspar_address}%s -O {output} --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"' % motifs[wildcards.condition])

rule motif_enrichment:
    input:
        fasta="results/{caller}/peak_fasta/{celltype}_{condition}.fa",
        meme="results/motifs/{condition}.meme"
    output:
        multiext("results/{caller}/motif_matches/{celltype}_{condition}/fimo", ".html", ".xml", ".tsv", ".gff") 
    shell:
        "fimo --oc results/{wildcards.caller}/motif_matches/{wildcards.celltype}_{wildcards.condition}/ {input.meme} {input.fasta}"

rule motif_plot:
    input:
        matches="results/{caller}/motif_matches/{sample}/fimo.tsv",
        peaks="results/{caller}/sorted_peaks/{sample}.narrowPeak"
    output:
        report("results/{caller}/motif_plots/{sample}.png", category="Motif_plots"),
        "results/{caller}/motif_plots/{sample}.pkl"
    script:
        "../scripts/motifplot.py"

rule join_plot:
    input:
        expand("results/{caller}/motif_plots/{{sample}}.pkl", caller=config["callers"])
    output:
        report("results/plots/motif/{sample}.png", category="Motif plots")
    script:
        "../scripts/joinplots.py"

rule sort_peaks:
    input:
        "results/{caller}/peaks/{filename}"
    output:
        "results/{caller}/sorted_peaks/{filename}"
    shell:
        "sort -nr -k5 {input} > {output}"

rule get_peak_sequences:
    input:
        peaks="results/{caller}/peaks/{sample}.narrowPeak",
    output:
        "results/{caller}/peak_fasta/{sample}.fa"
    params:
        reference=config["reference_fasta"]
    shell:
        "bedtools getfasta -fi {params.reference} -bed {input.peaks} > {output}"

