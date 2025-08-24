
rule kraken_uniq:
    """
    Performs Kraken Uniq classification

    Note that for the majority of samples, the time can probably 
    be reduced to 120 minutes and 40G of memory. 
    Typically, use 5 threads for this rule. The program does not efficiently use threads.
    """
    input: 
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        r_silva = RESULTS + "MTX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq"
    output:
        kraken_out = RESULTS + "MTX/kraken_uniq/{sample}/{sample}_krakenuniq_out.tsv",
        kraken_report = RESULTS + "MTX/kraken_uniq/{sample}/{sample}_krakenuniq_report.tsv"
    threads: 15
    resources:
        mem="120G",
        runtime = 480,
        slurm_partition="genomics"
    params:
        kraken_db = RESULTS + "databases/krakenuniq_micro/"
    conda:
        "../envs/krakenuniq.yaml" 
    shell:
        """
        krakenuniq -db {params.kraken_db} --threads {threads} --report-file {output.kraken_report} --preload-size 15G --output {output.kraken_out} --fastq-input {input.r1} {input.r_silva}
        """ 

use rule bracken as bracken_uniq with:
    input: 
        kraken_report = RESULTS + "MTX/kraken_uniq/{sample}/{sample}_krakenuniq_report.tsv"
    output:
        bracken_out = RESULTS + "MTX/kraken_uniq/{sample}/{sample}.bracken.tsv",
        bracken_report = RESULTS + "MTX/kraken_uniq/{sample}/{sample}.breport.tsv"
    params:
        read_length = "75",
        taxonomic_level = "S",
        read_threshold = "50",
        kraken_db = RESULTS + "databases/krakenuniq_micro/"
    resources:
        mem="25G",
        runtime = 60

use rule merge_bracken as merge_bracken_uniq with:
    input: 
        bracken_out = expand(RESULTS + "MTX/kraken_uniq/{sample}/{sample}.bracken.tsv", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports,
    output:
        merged_reports = RESULTS + "MTX/kraken_uniq/merged_bracken_uniq_report_profile.tsv"



rule kraken_uniq_exact:
    """
    Performs Kraken Uniq classification

    Note that for the majority of samples, the time can probably 
    be reduced to 120 minutes and 40G of memory. 
    Typically, use 5 threads for this rule. The program does not efficiently use threads.
    """
    input: 
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        r_silva = RESULTS + "MTX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq"
    output:
        kraken_out = RESULTS + "MTX/kraken_uniq2/{sample}/{sample}_krakenuniq_out.tsv",
        kraken_report = RESULTS + "MTX/kraken_uniq2/{sample}/{sample}_krakenuniq_report.tsv"
    threads: 25
    resources:
        mem="100G",
        runtime = 480,
        slurm_partition="genomics"
    params:
        kraken_db = RESULTS + "databases/krakenuniq_micro/"
    conda:
        "../envs/krakenuniq.yaml" 
    shell:
        """
        krakenuniq -db {params.kraken_db} --threads {threads} --report-file {output.kraken_report} --preload-size 15G --output {output.kraken_out} --fastq-input {input.r1} {input.r_silva}
        """ 

use rule merge_kraken as merge_kraken_uniq with:
    input: 
        kraken_reports = expand(RESULTS + "MTX/kraken_uniq/{sample}/{sample}_krakenuniq_report.tsv", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports
    output:
        merged_reports = RESULTS + "MTX/kraken_uniq/merged_kraken_uniq_report_profile.tsv"
    params:
        sample_list = list(samples_MTX["sample"])

