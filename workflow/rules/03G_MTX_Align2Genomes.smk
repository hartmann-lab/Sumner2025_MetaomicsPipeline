rule align_contigs:
    input:
        votus="/projects/b1042/HartmannLab/lung_phage/derep_viruses/vOTUs_genomad.fna",
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq"
    output:
        bam = RESULTS + "MTX/contig_alignment/{sample}/{sample}.bam"
    threads: 20
    shell:
        """
        module load minimap2
        module load samtools
        mkdir -p $(dirname {output.bam})

        minimap2 -x sr -a -t 20 {input.votus} {input.r1} | \
        samtools view -bS -F4 - | samtools sort -@ 16 - -o {output.bam} 
        """


rule coverm:
    input:
        bam = RESULTS + "MTX/contig_alignment/{sample}/{sample}.bam"
    output:
        coverage = RESULTS + "MTX/contig_alignment/{sample}/{sample}.tsv"
    threads: 1
    resources:
        mem = "5G",
        time="00:20:00"
    conda:
        "../envs/coverm.yml"
    shell:
        """
        coverm contig -b {input.bam} \
        -m mean rpkm tpm covered_fraction \
        --min-read-percent-identity 0.9 --min-covered-fraction 0.4 \
        -o {output.coverage}
        """
