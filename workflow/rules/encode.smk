rule encode_download:
    output:
        "results/reads/{sample}.bam"
    shell:
        temp("wget https://www.encodeproject.org/files/{wildcards.sample}/@@download/{wildcards.sample}.bam -O {output}")

def region_filter():
    if "chrom" not in config:
        return "1"
    f = '$1=="%s"' % config["chrom"]
    if "max_pos" in config:
        f += " && $3<%s" % config["max_pos"]
    return f
        
region_f = region_filter()
print(region_f)
rule region_bamtobed:
    input:
        "results/reads/{sample}.bam",
    output:
        temp("results/unsorted_reads/{sample}.bed")
    shell:
        """bedtools bamtobed -i {input} | awk '{{if ({region_f}) print}}' > {output}"""


rule sortbed:
    input:
        "results/unsorted_reads/{sample}.bed"
    output:
        "results/reads/{sample}.bed.gz"
    shell:
        "sort -k 1,1 -k2,2n {input} | gzip > {output}"
