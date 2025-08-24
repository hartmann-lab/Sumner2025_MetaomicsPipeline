
rule get_bowtie_human_idx:
    output:
        directory = directory(RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY/"),
        index = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY/chm13.draft_v1.0_plusY.1.bt2",
    shell:
        """
        cd $(dirname {output.directory})
        wget https://genome-idx.s3.amazonaws.com/bt/chm13.draft_v1.0_plusY.zip
        unzip chm13.draft_v1.0_plusY.zip
        """


rule get_fasta_index_bowtie2:
    input:
        index = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY/chm13.draft_v1.0_plusY.1.bt2"
    output:
        fasta = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY.fasta",
        faidx = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY.fasta.fai"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        bowtie2-inspect $(dirname {input.index})/$(basename {input.index} .1.bt2) > {output.fasta}
        samtools faidx {output.fasta}
        """


rule get_bowtie_silva_idx:
    output:
        directory = directory(RESULTS + "databases/bowtie2_silva/"),
        index = RESULTS + "databases/bowtie2_silva/SILVA_128_LSUParc_SSUParc_ribosomal_RNA.1.bt2l"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        kneaddata_database --download ribosomal_RNA bowtie2 {output.directory}
        """


rule fix_reads_spaces:
    """
    For using kneaddata with headers with '@VH00481:160:AAFLVL2M5:1:1101:65532:1114 1:N:0:GCATGT' 
    sourced one liner from biobakery forum post https://forum.biobakery.org/t/all-paired-end-read-unmatched/2895/35
    """
    input:
        r1 = get_mgx_read1,
        r2 = get_mgx_read2
    output:
        r1 = RESULTS + "MGX/reads/{sample}_R1.fastq",
        r2 = RESULTS + "MGX/reads/{sample}_R2.fastq"
    shell:
        """
        if [[ {input.r1} == *.fastq.gz ]]
        then
            echo "file is gzipped"
            zcat {input.r1} | awk 'NR%4==1 {{gsub(/ 1/, ""); print $0"#0/1"; next}} {{print}}' > {output.r1}
            zcat {input.r2} | awk 'NR%4==1 {{gsub(/ 2/, ""); print $0"#0/2"; next}} {{print}}' > {output.r2}
        elif [[ {input.r1} == *.fastq ]]
        then
            echo "file is not gzipped"
            cat {input.r1} | awk 'NR%4==1 {{gsub(/ 1/, ""); print $0"#0/1"; next}} {{print}}' > {output.r1}
            cat {input.r2} | awk 'NR%4==1 {{gsub(/ 2/, ""); print $0"#0/2"; next}} {{print}}' > {output.r2}
        else
            echo "file is not nothing"
        fi
        """


rule kneaddata:
    input:
        r1 = rules.fix_reads_spaces.output.r1,
        r2 = rules.fix_reads_spaces.output.r2,
        index_human = rules.get_bowtie_human_idx.output.index,
        index_silva = rules.get_bowtie_silva_idx.output.index

    output:
        RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq",
        RESULTS + "MGX/kneaddata/{sample}/{sample}_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_paired_contam_1.fastq",
        RESULTS + "MGX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.sam",
        #RESULTS + "MGX/kneaddata/{sample}_paired_1.fastq",
        #RESULTS + "MGX/kneaddata/{sample}_paired_2.fastq",
        #RESULTS + "MGX/kneaddata/{sample}_unmatched_1.fastq",
        #RESULTS + "MGX/kneaddata/{sample}_unmatched_2.fastq"
    threads: 5
    resources:
        mem="50G",
        runtime=220
    log: RESULTS + "MGX/{sample}.log"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        JAR_DIR=$(dirname $(which kneaddata))/../share/trimmomatic-0.33-3
        kneaddata --input1 {input.r1} --input2 {input.r2} --output $(dirname {output[0]}) --threads {threads} --output-prefix {wildcards.sample} --reference-db  $(dirname {input.index_human}) --reference-db  $(dirname {input.index_silva}) --sequencer-source TruSeq3 --run-fastqc-start --run-fastqc-end --cat-final-output --run-trim-repetitive --serial --store-temp-output
        """
#--remove-intermediate-output

rule kneaddata_count_table:
    input:
        expand(RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq", sample=samples_MGX["sample"])
    output:
        RESULTS + "MGX/kneaddata/mgx_kneaddata_count_table.tsv"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        mkdir $(dirname {output[0]})/logs

        echo "{input}" | sed "s/.fastq/.log/g" | xargs cp -t $(dirname {output})/logs

        kneaddata_read_count_table --input $(dirname {output[0]})/logs --output {output}

        rm $(dirname {output[0]})/logs/*
        rmdir $(dirname {output[0]})/logs
        """


rule multiqc:
    input:
        expand(RESULTS + "MGX/kneaddata/{sample}/{sample}.fastq", sample=samples_MGX["sample"])
    output:
        RESULTS + "MGX/multiqc/mgx_kneaddata_fastqc_report.html"
    params:
        search_pattern = RESULTS + "MGX/kneaddata/*/fastqc"
    conda: "../envs/multiqc.yaml"
    shell:
        """
        multiqc {params.search_pattern} \
        --filename $(basename {output[0]} .html0 \
        -o $(dirname {output[0]}) \
        --force
        """



rule csome_count:
    """
    Count the number of reads mapped to the each chromosome mapped chm13.draft_v1.0_plusY genome
    Takes SAM file generated by kneaddata, and makes sorted bam file for counting, and then counts the number of reads mapped to each chromosome. Finally, adds sample name to the beginning of each line for easier merging.
    """
    input:
        sam = RESULTS + "MGX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.sam",
        faidx = RESULTS + "databases/bowtie2/chm13.draft_v1.0_plusY.fasta.fai"
    output:
        bam = RESULTS + "MGX/kneaddata/{sample}/{sample}_chm13.draft_v1.0_plusY_bowtie2.bam",
        table = RESULTS + "MGX/kneaddata/{sample}/{sample}_chromosomal_coverage.tsv"
    conda: "../envs/kneaddata.yaml"
    shell:
        """
        samtools view -bS {input.sam} -t {input.faidx} | samtools sort -o {output.bam}
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.table}
        sed -i -e 's/^/{wildcards.sample}\t/' {output.table}
        """


rule merge_chromosomal_coverage:
    input:
        expand(RESULTS + "MGX/kneaddata/{sample}/{sample}_chromosomal_coverage.tsv", sample=samples_MGX["sample"])
    output:
        RESULTS + "MGX/kneaddata/mgx_human_chromosomal_coverage.tsv"
    shell:
        """
        cat {input} > {output[0]}
        sed -i -e 1i"sample\tfeature\tfeature_length\tread_count\tunmapped" {output[0]}
        sed -i -e "s/\\*/unplaced/" {output[0]}
        """
