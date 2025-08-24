rule get_bowtie_transctiptome_idx:
    output:
        directory = directory(RESULTS + "databases/bowtie2_transctiptome/"),
        index = RESULTS + "databases/bowtie2_transctiptome/human_hg38_refMrna.1.bt2"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        kneaddata_database --download human_transcriptome bowtie2 {output.directory}
        """


rule fix_reads_spaces_mtx:
    """
    For using kneaddata with headers with '@VH00481:160:AAFLVL2M5:1:1101:65532:1114 1:N:0:GCATGT' 
    sourced one liner from biobakery forum post https://forum.biobakery.org/t/all-paired-end-read-unmatched/2895/35
    """
    input:
        r1 = get_mtx_read1
    output:
        r1 = RESULTS + "MTX/reads/{sample}_R1.fastq",
    shell:
        """
        if [[ {input.r1} == *.fastq.gz ]]
        then
            echo "file is gzipped"
            zcat {input.r1} | awk 'NR%4==1 {{gsub(/ 1/, ""); print $0"#0/1"; next}} {{print}}' > {output.r1}
        elif [[ {input.r1} == *.fastq ]]
        then
            echo "file is not gzipped"
            cat {input.r1} | awk 'NR%4==1 {{gsub(/ 1/, ""); print $0"#0/1"; next}} {{print}}' > {output.r1}
        else
            echo "file is not nothing"
        fi
        """


rule kneaddata_se:
    input:
        r1 = rules.fix_reads_spaces_mtx.output.r1,
        index_human = rules.get_bowtie_human_idx.output.index,
        index_silva = rules.get_bowtie_silva_idx.output.index,
        index_transcriptome = rules.get_bowtie_transctiptome_idx.output.index

    output:
        RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq",
        RESULTS + "MTX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq",
        RESULTS + "MTX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.sam"
    threads: 15
    params:
        adapters = RESULTS + "MTX/kneaddata/{sample}/adapters.fa"
    resources:
        mem="20g",
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        JAR_DIR=$(dirname $(which kneaddata))/../share/trimmomatic-0.33-3
        kneaddata --unpaired {input.r1} --output $(dirname {output[0]}) --threads {threads} --output-prefix {wildcards.sample} --reference-db  $(dirname {input.index_human}) --reference-db  $(dirname {input.index_transcriptome}) --reference-db  $(dirname {input.index_silva}) --sequencer-source TruSeq3 --run-fastqc-start --run-trim-repetitive --run-fastqc-end --serial --store-temp-output
        """
        #kneaddata --unpaired {input.r1} --output $(dirname {output[0]}) --threads {threads} --output-prefix {wildcards.sample} --reference-db  $(dirname {input.index_human}) --reference-db  $(dirname {input.index_silva}) --sequencer-source TruSeq3 --run-fastqc-start --run-fastqc-end  --run-trim-repetitive --serial --store-temp-output --trimmomatic-options="ILLUMINACLIP:{params.adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 HEADCROP:10"


use rule kneaddata_count_table as kneaddata_count_table_se with:
    input:
        expand(RESULTS + "MTX/kneaddata/{sample}/{sample}.fastq", sample=samples_MTX["sample"])
    output:
        RESULTS + "MTX/kneaddata/mtx_kneaddata_count_table.tsv"
    conda: "../envs/kneaddata.yaml"


use rule csome_count as csome_count_mtx with:
    input:
        sam = RESULTS + "MTX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.sam",
        faidx = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY.fasta.fai"
    output:
        bam = RESULTS + "MTX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.bam",
        table = RESULTS + "MTX/kneaddata/{sample}/{sample}_chromosomal_coverage.tsv"


use rule merge_chromosomal_coverage as merge_chromosomal_coverage_mtx with:
    input:
        expand(RESULTS + "MTX/kneaddata/{sample}/{sample}_chromosomal_coverage.tsv", sample=samples_MTX["sample"])
    output:
        RESULTS + "MTX/kneaddata/mtx_human_chromosomal_coverage.tsv"
