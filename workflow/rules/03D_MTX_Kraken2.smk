
rule get_kraken2_db:
    output:
        file = RESULTS + "databases/kraken_silva/16S_Silva138_20200326.tgz",
        kraken2_db =  directory(RESULTS + "databases/kraken_silva/16S_SILVA138_k2db")
    params:
        kraken2_db = config["KRAKEN_DATABASE_SILVA"]
    shell:
        """
        cd $(dirname {output.file})
        wget {params.kraken2_db}
        tar -xvzf 16S_Silva138_20200326.tgz -C ./
        """


rule get_kraken_tools:
    output:
        file = RESULTS + "databases/kraken_tools/" + config["KRAKEN_TOOLS"][1] + ".zip",
        kraken_tools = directory(RESULTS + "databases/kraken_tools/KrakenTools-" + config["KRAKEN_TOOLS"][1]),
        script_combine_kreports = RESULTS + "databases/kraken_tools/KrakenTools-" + config["KRAKEN_TOOLS"][1] + "/combine_kreports.py",
        script_kreport2mpa = RESULTS + "databases/kraken_tools/KrakenTools-" + config["KRAKEN_TOOLS"][1] + "/kreport2mpa.py"
    params:
        url = config["KRAKEN_TOOLS"][0],
        commit = config["KRAKEN_TOOLS"][1]
    shell:
        """
        cd $(dirname {output.file})
        wget {params.url}
        unzip $(basename {output.file})
        """


rule kraken2_silva_mtx:
    """
    Performs taxnomic classification with Kraken2 

    Outputs a kraken2-style report and metaphlan-style report with a
    script from KrakenTools
    """
    input: 
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq",
    output:
        kraken_out = RESULTS + "MTX/kraken2_silva/{sample}/{sample}_kraken2out.txt",
        kraken_report = RESULTS + "MTX/kraken2_silva/{sample}/{sample}_kraken2report.txt"
    threads: 19
    resources:
        mem="60G",
        runtime = 45
    params:
        kraken_db = RESULTS + "databases/kraken_silva/16S_SILVA138_k2db", #/software/kraken/database/kraken_db
        minimum_hit_groups = 3,
        extra_reads =  ""
    shell:
        """
        module load kraken/2
        kraken2 --threads {threads} \
            --output {output.kraken_out} \
            --report {output.kraken_report} \
            --db {params.kraken_db} \
            --minimum-hit-groups {params.minimum_hit_groups} \
            {input.r1} {params.extra_reads} 
        """ 


rule merge_kraken:
    """
    Outputs a metaphlan-style report with a script from KrakenTools
    """
    input: 
        kraken_reports = expand(RESULTS + "MTX/kraken2_silva/{sample}/{sample}_kraken2report.txt", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports
    output:
        merged_reports = RESULTS + "MTX/kraken2_silva/merged_kraken_report_profile.tsv"
    params:
        sample_list = list(samples_MTX["sample"])
    threads: 1
    resources:
        mem="15G",
        runtime = 60
    shell:
        """
        {input.script} -r {input.kraken_reports} -o {output.merged_reports} --sample-names {params.sample_list} --display-headers 
        """ 


rule bracken:
    """
    Performs abundance estimation from with Kraken2 classification
    """
    input: 
        kraken_report = RESULTS + "MTX/kraken2_silva/{sample}/{sample}_kraken2report.txt"
    output:
        bracken_out = RESULTS + "MTX/kraken2_silva/{sample}/{sample}.bracken.tsv",
        bracken_report = RESULTS + "MTX/kraken2_silva/{sample}/{sample}.breport.tsv"
    threads: 1
    conda:
        "../envs/bracken.yaml"
    resources:
        mem="5G",
        runtime = 10
    params:
        read_length = "75",
        taxonomic_level = "G",
        read_threshold = "100",
        kraken_db = RESULTS + "databases/kraken_silva/16S_SILVA138_k2db" #/software/kraken/database/kraken_db
    shell:
        """
        module load kraken/2
        bracken -d {params.kraken_db} \
            -i {input.kraken_report} \
            -r {params.read_length} \
            -l {params.taxonomic_level} \
            -t {params.read_threshold} \
            -o {output.bracken_out} \
            -w {output.bracken_report} \
        """


rule merge_bracken:
    input: 
        bracken_out = expand(RESULTS + "MTX/kraken2_silva/{sample}/{sample}.bracken.tsv", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports,
    output:
        merged_reports = RESULTS + "MTX/kraken2_silva/merged_bracken_report_profile.tsv"
    conda:
        "../envs/bracken.yaml"
    resources:
        mem="15G",
        runtime = 30
    shell:
        """
        combine_bracken_outputs.py --files {input.bracken_out} --output {output.merged_reports}
        """


# STANDARD DATABASE
"""
Following code is implementation of above with standard database instead of silva
"""

use rule kraken2_silva_mtx as kraken2_standard_mtx with:
    input: 
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
    output:
        kraken_out = RESULTS + "MTX/kraken2_standard/{sample}/{sample}_kraken2out.txt",
        kraken_report = RESULTS + "MTX/kraken2_standard/{sample}/{sample}_kraken2report.txt"
    threads: 19
    resources:
        mem="100G",
        runtime = 15
    params:
        kraken_db = RESULTS + "databases/kraken_standard", #/software/kraken/database/kraken_db
        minimum_hit_groups = 3,
        extra_reads =  RESULTS + "MTX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq"



use rule merge_kraken as merge_kraken_standard with:
    input: 
        kraken_reports = expand(RESULTS + "MTX/kraken2_standard/{sample}/{sample}_kraken2report.txt", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports
    output:
        merged_reports = RESULTS + "MTX/kraken2_standard/merged_kraken_report_profile.tsv"
    params:
        sample_list = list(samples_MTX["sample"])


use rule bracken as bracken_standard with:
    input: 
        kraken_report = RESULTS + "MTX/kraken2_standard/{sample}/{sample}_kraken2report.txt"
    output:
        bracken_out = RESULTS + "MTX/kraken2_standard/{sample}/{sample}.bracken.tsv",
        bracken_report = RESULTS + "MTX/kraken2_standard/{sample}/{sample}.breport.tsv"
    params:
        read_length = "75",
        taxonomic_level = "G",
        read_threshold = "100",
        kraken_db = RESULTS + "databases/kraken_standard" #/software/kraken/database/kraken_db

use rule merge_bracken as merge_bracken_standard with:
    input: 
        bracken_out = expand(RESULTS + "MTX/kraken2_standard/{sample}/{sample}.bracken.tsv", sample=samples_MTX["sample"]),
        file = rules.get_kraken_tools.output.kraken_tools,
        script = rules.get_kraken_tools.output.script_combine_kreports,
    output:
        merged_reports = RESULTS + "MTX/kraken2_standard/merged_bracken_report_profile.tsv"


