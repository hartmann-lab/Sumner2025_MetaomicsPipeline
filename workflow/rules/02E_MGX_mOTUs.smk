rule get_motus_idx:
    output:
        index = RESULTS + "databases/motus/motus_index"
    threads: 20
    conda:
        "../envs/motus.yml"
    shell:
        """
        motus downloadDB 
        mkdir -p $(dirname {output.index})
        touch {output.index}
        """


rule motus:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq",
        index = rules.get_motus_idx.output.index
    output:
        RESULTS + "MGX/motus/{sample}/{sample}.tsv"
    params:
        min_alingnment_length = 75
    threads: 20
    resources:
        mem = "30G"
    conda:
        "../envs/motus.yml"
    shell:
        """
        motus profile -s {input.r1} -n {wildcards.sample} -o {output} -t {threads} -l {params.min_alingnment_length}
        """


rule motus_count_table:
    input:
        expand(RESULTS + "MGX/motus/{sample}/{sample}.tsv", sample=samples_MGX["sample"])
    output:
        RESULTS + "MGX/motus/mgx_motus_count_table.tsv"
    params:
        input_list = ",".join(expand(RESULTS + "MGX/motus/{sample}/{sample}.tsv", sample=samples_MGX["sample"]))
    conda:
        "../envs/motus.yml"
    shell:
        """
        motus merge {params.input_list} -o {output} -a human,air,freshwater,mouse
        """ 

