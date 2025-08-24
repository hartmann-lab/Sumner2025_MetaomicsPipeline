use rule metaphlan as metaphlan_se with:
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db_file,
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq"
    output:
        profile =  RESULTS + "MTX/metaphlan/{sample}/{sample}_profile.txt",
        vsc_profile =  RESULTS + "MTX/metaphlan/{sample}/{sample}_vsc.txt",
        bowtie_out =  RESULTS + "MTX/metaphlan/{sample}/{sample}.bowtie2.bz2"
    params:
        metaphlan_idx = config["METAPHLAN_DATABASE"][0], # Index for metaphlan
        read_min_len = 45,
        stat_q = 0.05



use rule metaphlan_merge as metaphlan_merge_mtx with:
    input:
        expand(RESULTS + "MTX/metaphlan/{sample}/{sample}_profile.txt", sample=samples_MTX["sample"])
    output:
        RESULTS + "MTX/metaphlan/mtx_metaphlan_abundance_table_all.txt"


use rule metaphlan_species_abundance as metaphlan_species_abundance_mtx with:
    input:
        RESULTS + "MTX/metaphlan/mtx_metaphlan_abundance_table_all.txt"
    output:
        species = RESULTS + "MTX/metaphlan/mtx_metaphlan_abundance_table_species.txt",
        genus = RESULTS + "MTX/metaphlan/mtx_metaphlan_abundance_table_genus.txt"
