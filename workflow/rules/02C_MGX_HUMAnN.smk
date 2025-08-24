rule get_humann_databases:
    output:
        directory = directory(RESULTS + "databases/humann")
    conda: 
        "../envs/humann.yaml"
    shell:
        """
        humann_databases --download chocophlan full {output.directory}
        humann_databases --download uniref uniref90_diamond {output.directory}
        """


rule get_humann_mapping:
    """
    Realized  needed this after installling above humann rule. 
    In future, combine this with the above rule.
    """
    output:
        mapping_directory = directory(RESULTS + "databases/humann_mapping")
    conda: 
        "../envs/humann.yaml"
    shell:
        """
        humann_databases --download utility_mapping full {output.mapping_directory}
        """ 


rule humann:
    input:
        database = directory(RESULTS + "databases/humann"),
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq",
        metaphlan_profile = RESULTS + "MGX/metaphlan/{sample}/{sample}_profile.txt"
    output:
        gene_fam = RESULTS + "MGX/humann/{sample}/{sample}_genefamilies.tsv",
        path_cov = RESULTS + "MGX/humann/{sample}/{sample}_pathcoverage.tsv",
        path_abund = RESULTS + "MGX/humann/{sample}/{sample}_pathabundance.tsv"
    params:
        translated_coverage_threshold = 50.0,
        extra_params = ""
    conda: 
        "../envs/humann.yaml"
    threads: 25
    resources:
        mem="70G",
        runtime=240
    shell:
        """
        humann --input {input.r1} \
        --output $(dirname {output.gene_fam}) \
        --threads {threads} \
        --taxonomic-profile {input.metaphlan_profile} \
        --translated-subject-coverage-threshold {params.translated_coverage_threshold} \
        --nucleotide-database {input.database}/chocophlan \
        --protein-database {input.database}/uniref \
        {params.extra_params}
        """


rule humann_merge:
    input:
        gf = expand(RESULTS + "MGX/humann/{sample}/{sample}_genefamilies.tsv", sample=samples_MGX["sample"]),
        pc = expand(RESULTS + "MGX/humann/{sample}/{sample}_pathcoverage.tsv", sample=samples_MGX["sample"]),
        pa = expand(RESULTS + "MGX/humann/{sample}/{sample}_pathabundance.tsv", sample=samples_MGX["sample"])
    output:
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies.tsv",
        pc = RESULTS + "MGX/humann/mgx_humann_pathcoverage.tsv",
        pa = RESULTS + "MGX/humann/mgx_humann_pathabundance.tsv"
    params:
        humann_dir = RESULTS + "MGX/humann/"
    conda: 
        "../envs/humann.yaml"
    shell:
        """
        humann_join_tables --input {params.humann_dir} --output {output.gf} --file_name genefamilies --search-subdirectories
        humann_join_tables --input {params.humann_dir} --output {output.pc} --file_name pathcoverage --search-subdirectories
        humann_join_tables --input {params.humann_dir} --output {output.pa} --file_name pathabundance --search-subdirectories
        """


rule humann_renorm_table:
    input:
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies.tsv",
        pc = RESULTS + "MGX/humann/mgx_humann_pathcoverage.tsv",
        pa = RESULTS + "MGX/humann/mgx_humann_pathabundance.tsv",
    output:
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm.tsv",
        pc = RESULTS + "MGX/humann/mgx_humann_pathcoverage_cpm.tsv",
        pa = RESULTS + "MGX/humann/mgx_humann_pathabundance_cpm.tsv"
    conda: "../envs/humann.yaml"
    resources:
        mem="120G",
        runtime=240
    shell:
        """
        humann_renorm_table --input {input.gf} --output {output.gf} --units relab --update-snames
        humann_renorm_table --input {input.pc} --output {output.pc} --units relab --update-snames
        humann_renorm_table --input {input.pa} --output {output.pa} --units relab --update-snames
        """


rule humann_regroup_table:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm_KEGGOrthology.tsv"
    conda: "../envs/humann.yaml"
    resources:
        mem="60G",
        runtime=240
    shell:
        """
        humann_regroup_table --input {input.gf} --output {output.gf} --groups uniref90_ko
        """


rule humann_regroup_table_pfam:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm_PFAM.tsv"
    conda: "../envs/humann.yaml"
    resources:
        mem="30G",
        runtime=480
    shell:
        """
        humann_regroup_table --input {input.gf} --output {output.gf} --groups uniref90_pfam
        """


rule humann_split_stratified_table:
    input:
        ko = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm_KEGGOrthology.tsv",
        pf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm_PFAM.tsv",
        gf = RESULTS + "MGX/humann/mgx_humann_genefamilies_cpm.tsv",
        pc = RESULTS + "MGX/humann/mgx_humann_pathcoverage_cpm.tsv",
        pa = RESULTS + "MGX/humann/mgx_humann_pathabundance_cpm.tsv"
    output:
        ko = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_PFAM_unstratified.tsv",
        gf = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_unstratified.tsv",
        pc = RESULTS + "MGX/final_tables/mgx_humann_pathcoverage_cpm_unstratified.tsv",
        pa = RESULTS + "MGX/final_tables/mgx_humann_pathabundance_cpm_unstratified.tsv"
    conda: "../envs/humann.yaml"
    params:
        humann_dir = RESULTS + "MGX/final_tables/"
    shell:
        """
        humann_split_stratified_table --input {input.ko} --output {params.humann_dir}
        humann_split_stratified_table --input {input.gf} --output {params.humann_dir}
        humann_split_stratified_table --input {input.pc} --output {params.humann_dir}
        humann_split_stratified_table --input {input.pa} --output {params.humann_dir}
        humann_split_stratified_table --input {input.pf} --output {params.humann_dir}
        """


rule humann_rename_table:
    input:
        ko = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_PFAM_unstratified.tsv"
    output:
        ko = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv",
        pf = RESULTS + "MGX/final_tables/mgx_humann_genefamilies_cpm_PFAM_unstratified_named.tsv",
    conda: "../envs/humann.yaml"
    resources:
        mem="30G",
        runtime=60
    shell:
        """
        humann_rename_table --input {input.ko} --output {output.ko} --names kegg-orthology
        humann_rename_table --input {input.pf} --output {output.pf} --names pfam
        """
