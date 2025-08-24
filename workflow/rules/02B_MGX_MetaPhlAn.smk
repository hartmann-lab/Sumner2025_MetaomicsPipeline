
rule metaphlan_setup:
    output:
        metaphlan_db=directory(RESULTS + "databases/chocophlan"),
        metaphlan_db_file = RESULTS + "databases/chocophlan/{}.rev.2.bt2l".format(config["METAPHLAN_DATABASE"][0])
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["METAPHLAN_DATABASE"][0] # Index for metaphlan
    threads: 5
    resources:
        mem="200g",
        time="04:00:00"
    shell:
        """
        metaphlan --install --index {params.metaphlan_idx} --bowtie2db {output.metaphlan_db} --nproc {threads}
        """


rule metaphlan:
    """
    Min mapq value & stat_q value are set to 1 and 0.1 respectively to reduce the number of unclassified reads.
    https://forum.biobakery.org/t/unexpected-difference-with-very-sensitive-local-for-environmental-sample/4305
    https://forum.biobakery.org/t/understanding-parameters-stat-q-for-environmental-sample/2204
    """
    input:
        metaphlan_db = rules.metaphlan_setup.output.metaphlan_db_file,
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq"
    output:
        profile =  RESULTS + "MGX/metaphlan/{sample}/{sample}_profile.txt",
        vsc_profile =  RESULTS + "MGX/metaphlan/{sample}/{sample}_vsc.txt",
        bowtie_out =  RESULTS + "MGX/metaphlan/{sample}/{sample}.bowtie2.bz2"
    conda: 
        "../envs/metaphlan.yml"
    params:
        metaphlan_idx = config["METAPHLAN_DATABASE"][0], # Index for metaphlan
        read_min_len = 70,
        stat_q = 0.1
    threads: 5
    resources:
        mem="30g",
        time="04:00:00"
    shell:
        """
        metaphlan {input.r1} \
        --bowtie2out {output.bowtie_out} \
        --index {params.metaphlan_idx} \
        --bowtie2db $(dirname {input.metaphlan_db}) \
        --nproc {threads} \
        --input_type fastq \
        --unclassified_estimation \
        -t rel_ab \
        -o {output.profile} \
        --profile_vsc \
        --vsc_out {output.vsc_profile} \
        --stat_q {params.stat_q} \
        --min_mapq_val 5 \
        --read_min_len {params.read_min_len}
        """

# --bt2_ps 'very-sensitive-local'  maybe? no this doesnt work


rule metaphlan_merge:
    input:
        expand(RESULTS + "MGX/metaphlan/{sample}/{sample}_profile.txt", sample=samples_MGX["sample"])
    output:
        RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
    conda: 
        "../envs/metaphlan.yml"
    shell:
        """
        merge_metaphlan_tables.py {input} > {output}
        """


rule metaphlan_species_abundance:
    input:
        RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_all.txt"
    output:
        species = RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_species.txt",
        genus = RESULTS + "MGX/metaphlan/mgx_metaphlan_abundance_table_genus.txt"

    shell:
        """
        echo "species level table"
        grep -E "t__|clade|UNCLASSIFIED" {input} | sed 's/^.*t__//g' | sed -e 's/clade_name/sample/g' > {output.species}

        echo "genus level table"
        grep -E "g__|clade|UNCLASSIFIED" {input} | grep -v "s__" | sed 's/^.*g__//g' | sed -e 's/clade_name/sample/g' > {output.genus}
        """
        
        
#grep -E "s__|clade|UNCLASSIFIED" {input} | grep -v "t__" | sed 's/^.*s__//g' | sed -e 's/clade_name/sample/g' > {output.species}
