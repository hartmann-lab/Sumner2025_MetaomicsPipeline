use rule humann as humann_se_mgxmpa_closed with:
    input:
        database = directory(RESULTS + "databases/humann"),
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        metaphlan_profile = RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
    output:
        gene_fam = RESULTS + "MTX/humann_closed/{sample}/{sample}_genefamilies.tsv",
        path_cov = RESULTS + "MTX/humann_closed/{sample}/{sample}_pathcoverage.tsv",
        path_abund = RESULTS + "MTX/humann_closed/{sample}/{sample}_pathabundance.tsv"
    params:
        translated_coverage_threshold = 25.0,
        extra_params = "--prescreen-threshold 0 --translated-identity-threshold 50 --translated-query-coverage-threshold 70 --nucleotide-query-coverage-threshold 70 --nucleotide-subject-coverage-threshold 25"



use rule humann_merge as humann_merge_se_mgxmpa_closed with:
    input:
        gf = expand(RESULTS + "MTX/humann_closed/{sample}/{sample}_genefamilies.tsv", sample=samples_MTX["sample"]),
        pc = expand(RESULTS + "MTX/humann_closed/{sample}/{sample}_pathcoverage.tsv", sample=samples_MTX["sample"]),
        pa = expand(RESULTS + "MTX/humann_closed/{sample}/{sample}_pathabundance.tsv", sample=samples_MTX["sample"])
    output:
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies.tsv",
        pc = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathcoverage.tsv",
        pa = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathabundance.tsv"
    params:
        humann_dir = RESULTS + "MTX/humann_closed/"


use rule humann_renorm_table as humann_renorm_table_se_mgxmpa_closed with:
    input:
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies.tsv",
        pc = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathcoverage.tsv",
        pa = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathabundance.tsv",
    output:
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm.tsv",
        pc = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathcoverage_cpm.tsv",
        pa = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathabundance_cpm.tsv"


use rule humann_regroup_table as humann_regroup_table_se_mgxmpa_closed with:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology.tsv"


use rule humann_regroup_table_pfam as humann_regroup_table_pfam_se_mgxmpa_closed with:
    input:
        mapping_directory = directory(RESULTS + "databases/humann_mapping"),
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm.tsv"
    output:
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM.tsv"


use rule humann_split_stratified_table as humann_split_stratified_table_se_mgxmpa_closed with:
    input:
        ko = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology.tsv",
        pf =RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM.tsv",
        gf = RESULTS + "MTX/humann_closed/mtx_humann_closed_genefamilies_cpm.tsv",
        pc = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathcoverage_cpm.tsv",
        pa = RESULTS + "MTX/humann_closed/mtx_humann_closed_pathabundance_cpm.tsv"
    output:
        ko = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM_unstratified.tsv",
        gf = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_unstratified.tsv",
        pc = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_pathcoverage_cpm_unstratified.tsv",
        pa = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_pathabundance_cpm_unstratified.tsv"
    params:
        humann_dir = RESULTS + "MTX/final_tables/humann_closed/"

use rule humann_rename_table as humann_rename_table_se_mgxmpa_closed with:
    input:
        ko = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology_unstratified.tsv",
        pf = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM_unstratified.tsv"
    output:
        ko = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_KEGGOrthology_unstratified_named.tsv",
        pf = RESULTS + "MTX/final_tables/humann_closed/mtx_humann_closed_genefamilies_cpm_PFAM_unstratified_named.tsv"




use rule humann as humann_se_mgxmpa_closed_rc with:
    input:
        database = directory(RESULTS + "databases/humann"),
        r1 = RESULTS + "rev_comp/{sample}_rc.fastq",
        metaphlan_profile = RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
    output:
        gene_fam = RESULTS + "MTX/humann_rc/{sample}/{sample}_genefamilies.tsv",
        path_cov = RESULTS + "MTX/humann_rc/{sample}/{sample}_pathcoverage.tsv",
        path_abund = RESULTS + "MTX/humann_rc/{sample}/{sample}_pathabundance.tsv"
    params:
        translated_coverage_threshold = 25.0,
        extra_params = "--prescreen-threshold 0 --translated-identity-threshold 50 --translated-query-coverage-threshold 70 --nucleotide-query-coverage-threshold 70 --nucleotide-subject-coverage-threshold 25 --diamond-options\"--top 1 --outfmt 6 --forwardonly\""


use rule humann as humann_se_mgxmpa_closed_chopped with:
    input:
        database = directory(RESULTS + "databases/humann"),
        r1 = RESULTS + "rev_comp/{sample}_chopped.fastq",
        metaphlan_profile = RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
    output:
        gene_fam = RESULTS + "MTX/humann_chopped/{sample}/{sample}_genefamilies.tsv",
        path_cov = RESULTS + "MTX/humann_chopped/{sample}/{sample}_pathcoverage.tsv",
        path_abund = RESULTS + "MTX/humann_chopped/{sample}/{sample}_pathabundance.tsv"
    params:
        translated_coverage_threshold = 50.0,
        extra_params = "--prescreen-threshold 0"
