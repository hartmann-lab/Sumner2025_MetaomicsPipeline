rule megahit:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        scaffolds = RESULTS + "MGX/megahit/{sample}/{sample}.contigs.fa",
        out_dir = directory(RESULTS + "MGX/megahit/{sample}")
    params:
        out_dir = RESULTS + "MGX/megahit/{sample}"
    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="30G",
        runtime=60
    shell:
        """
        megahit -t {threads} \
            -m 0.9 -r {input.r1} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """

rule megahit_paired:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        scaffolds = RESULTS + "MGX/megahit_paired/{sample}/{sample}.contigs.fa",
        out_dir = directory(RESULTS + "MGX/megahit_paired/{sample}")
    params:
        out_dir = RESULTS + "MGX/megahit_paired/{sample}",
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_1.fastq", 
        r2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_2.fastq", 
        r1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_1.fastq", 
        r2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_2.fastq", 
    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="30G",
        runtime=60
    shell:
        """
        megahit -t {threads} \
            -m 0.9 \
            -1 {params.r1} -2 {params.r2} \
            -r {params.r1_u},{params.r2_u} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """


rule megahit_paired2:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        scaffolds = RESULTS + "MGX/megahit_paired2/{sample}/{sample}.contigs.fa",
        out_dir = directory(RESULTS + "MGX/megahit_paired2/{sample}")
    params:
        out_dir = RESULTS + "MGX/megahit_paired2/{sample}",
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_1.fastq", 
        r2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_2.fastq", 
        r1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_1.fastq", 
        r2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_2.fastq", 
        
        rr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.1.fastq", 
        rr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.2.fastq", 
        rr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.1.fastq", 
        rr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.2.fastq", 
        
        sr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq", 
        sr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq", 
        sr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq", 
        sr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq"

    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="30G",
        runtime=360
    shell:
        """
        megahit -t {threads} \
            -m 0.9 \
            -1 {params.r1},{params.rr1},{params.sr1} -2 {params.r2},{params.rr2},{params.sr2} \
            -r {params.r1_u},{params.r2_u},{params.rr1_u},{params.rr2_u},{params.sr1_u},{params.sr2_u} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """


# 16s 
rule megahit_paired3:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        scaffolds = RESULTS + "MGX/megahit_paired3/{sample}/{sample}.contigs.fa",
        out_dir = directory(RESULTS + "MGX/megahit_paired3/{sample}")
    params:
        out_dir = RESULTS + "MGX/megahit_paired3/{sample}",
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_1.fastq", 
        r2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_2.fastq", 
        r1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_1.fastq", 
        r2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_2.fastq", 
        
        rr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.1.fastq", 
        rr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.2.fastq", 
        rr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.1.fastq", 
        rr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.2.fastq", 
        
        sr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq", 
        sr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq", 
        sr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq", 
        sr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq"

    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="30G",
        runtime=360
    shell:
        """
        megahit -t {threads} \
            -m 0.9 \
            -1 {params.r1},{params.sr1} -2 {params.r2},{params.sr2} \
            -r {params.r1_u},{params.r2_u},{params.sr1_u},{params.sr2_u} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """



# repetitive 
rule megahit_paired4:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        scaffolds = RESULTS + "MGX/megahit_paired4/{sample}/{sample}.contigs.fa",
        out_dir = directory(RESULTS + "MGX/megahit_paired4/{sample}")
    params:
        out_dir = RESULTS + "MGX/megahit_paired4/{sample}",
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_1.fastq", 
        r2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_paired_2.fastq", 
        r1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_1.fastq", 
        r2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_unmatched_2.fastq", 
        
        rr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.1.fastq", 
        rr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.2.fastq", 
        rr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.1.fastq", 
        rr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}.repeats.removed.unmatched.2.fastq", 
        
        sr1 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq", 
        sr2 = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_2.fastq", 
        sr1_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_1_contam.fastq", 
        sr2_u = RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_unmatched_2_contam.fastq"

    conda:
        "../envs/megahit.yml"
    threads: 20
    resources:
        mem="30G",
        runtime=360
    shell:
        """
        megahit -t {threads} \
            -m 0.9 \
            -1 {params.r1},{params.rr1} -2 {params.r2},{params.rr2} \
            -r {params.r1_u},{params.r2_u},{params.rr1_u},{params.rr2_u} \
            --out-prefix {wildcards.sample} \
            -o {params.out_dir}_tmp
        mv {params.out_dir}_tmp/* {params.out_dir}
        rmdir {params.out_dir}_tmp
        """