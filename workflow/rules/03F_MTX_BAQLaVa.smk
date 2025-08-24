
use rule baqlava as baqlava_mtx with:
    input:
        r1 = RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        metaphlan_profile = RESULTS + "MTX/metaphlan/{sample}/{sample}_profile.txt",
        index_prtn = rules.get_baqlava_idx.output.index_prtn,
        index_nuc = rules.get_baqlava_idx.output.index_nuc
    output:
        profile = RESULTS + "MTX/baqlava/{sample}/{sample}_BAQLaVa_profile.txt",
        out_dir = directory(RESULTS + "MTX/baqlava/{sample}/")
    params:
        extra="--bypass-bacterial-depletion"
    threads: 20
    resources:
        mem = "30G",
        time="01:00:00"
