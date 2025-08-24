use rule humann as humann_se with:
    input:
        database = directory(RESULTS + "databases/humann"),
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        metaphlan_profile = RESULTS + "MTX/metaphlan/{sample}/{sample}_profile.txt"
    output:
        gene_fam = RESULTS + "MTX/humann/{sample}/{sample}_genefamilies.tsv",
        path_cov = RESULTS + "MTX/humann/{sample}/{sample}_pathcoverage.tsv",
        path_abund = RESULTS + "MTX/humann/{sample}/{sample}_pathabundance.tsv"
    threads: 10
    resources:
        mem="50G",
        runtime=240
    params:
        translated_coverage_threshold = 25.0,
        extra_params = "--translated-identity-threshold 50 --translated-query-coverage-threshold 70 --nucleotide-query-coverage-threshold 70 --nucleotide-subject-coverage-threshold 25"



use rule humann_merge as humann_merge_se with:
    input:
        gf = expand(RESULTS + "MTX/humann/{sample}/{sample}_genefamilies.tsv", sample=samples_MTX["sample"]),
        pc = expand(RESULTS + "MTX/humann/{sample}/{sample}_pathcoverage.tsv", sample=samples_MTX["sample"]),
        pa = expand(RESULTS + "MTX/humann/{sample}/{sample}_pathabundance.tsv", sample=samples_MTX["sample"])
    output:
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies.tsv",
        pc = RESULTS + "MTX/humann/mtx_humann_pathcoverage.tsv",
        pa = RESULTS + "MTX/humann/mtx_humann_pathabundance.tsv"
    params:
        humann_dir = RESULTS + "MTX/humann/"


use rule humann_renorm_table as humann_renorm_table_se with:
    input:
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies.tsv",
        pc = RESULTS + "MTX/humann/mtx_humann_pathcoverage.tsv",
        pa = RESULTS + "MTX/humann/mtx_humann_pathabundance.tsv",
    output:
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm.tsv",
        pc = RESULTS + "MTX/humann/mtx_humann_pathcoverage_cpm.tsv",
        pa = RESULTS + "MTX/humann/mtx_humann_pathabundance_cpm.tsv"


use rule humann_regroup_table as humann_regroup_table_se with:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm_KEGGOrthology.tsv"

use rule humann_regroup_table_pfam as humann_regroup_table_pfam_se with:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm_PFAM.tsv"


use rule humann_split_stratified_table as humann_split_stratified_table_se with:
    input:
        ko = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm_KEGGOrthology.tsv",
        pf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm_PFAM.tsv",
        gf = RESULTS + "MTX/humann/mtx_humann_genefamilies_cpm.tsv",
        pc = RESULTS + "MTX/humann/mtx_humann_pathcoverage_cpm.tsv",
        pa = RESULTS + "MTX/humann/mtx_humann_pathabundance_cpm.tsv"
    output:
        ko = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_PFAM_unstratified.tsv",
        gf = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_unstratified.tsv",
        pc = RESULTS + "MTX/final_tables/mtx_humann_pathcoverage_cpm_unstratified.tsv",
        pa = RESULTS + "MTX/final_tables/mtx_humann_pathabundance_cpm_unstratified.tsv"
    params:
        humann_dir = RESULTS + "MTX/final_tables/"


use rule humann_rename_table as humann_rename_table_se with:
    input:
        ko = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_PFAM_unstratified.tsv"
    output:
        ko = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv",
        pf = RESULTS + "MTX/final_tables/mtx_humann_genefamilies_cpm_PFAM_unstratified_named.tsv"


