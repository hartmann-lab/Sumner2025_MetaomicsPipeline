
use rule motus as motus_mtx with:
    input:
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        index = rules.get_motus_idx.output.index
    output:
        RESULTS + "MTX/motus/{sample}/{sample}.tsv"
    params:
        min_alingnment_length = 40
    threads: 20


use rule motus_count_table as motus_count_table_mtx with:
    input:
        expand(RESULTS + "MTX/motus/{sample}/{sample}.tsv", sample=samples_MTX["sample"])
    output:
        RESULTS + "MTX/motus/mtx_motus_count_table.tsv"
    params:
        input_list = ",".join(expand(RESULTS + "MTX/motus/{sample}/{sample}.tsv", sample=samples_MTX["sample"]))