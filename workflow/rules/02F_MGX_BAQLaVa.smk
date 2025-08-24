rule get_baqlava_idx:
    output:
        index_prtn = RESULTS + "databases/baqlava/data/BAQLaVa.V0.5.protein/BAQLaVa_ORFs_90ID_201901b.dmnd",
        index_nuc = RESULTS + "databases/baqlava/data/BAQLaVa.V0.5.nucleotide/allmarkers_201901b.1.bt2",
        db = directory(RESULTS + "databases/baqlava")
    threads: 20
    conda:
        "../envs/baqlava.yml"
    shell:
        """
        mkdir -p $(dirname {output.index})
        baqlava_database --download database baqlava-db $(dirname {output.index})
        touch {output.index}
        """

rule baqlava:
    input:
        r1 = RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq",
        metaphlan_profile = RESULTS + "MGX/metaphlan/{sample}/{sample}_profile.txt",
        index_prtn = rules.get_baqlava_idx.output.index_prtn,
        index_nuc = rules.get_baqlava_idx.output.index_nuc
    output:
        profile = RESULTS + "MGX/baqlava/{sample}/{sample}_BAQLaVa_profile.txt",
        out_dir = directory(RESULTS + "MGX/baqlava/{sample}/")
    params:
        extra=""
    threads: 20
    resources:
        mem = "30G",
        time="01:00:00" 
    conda:
        "../envs/baqlava.yml"
    shell:
        """
        baqlava -i {input.r1} -o {output.out_dir} --threads {threads} --taxonomic-profile {input.metaphlan_profile} --protdb $(dirname {input.index_prtn}) --nucdb $(dirname {input.index_nuc}) --keep-tempfiles {params.extra}
        """

