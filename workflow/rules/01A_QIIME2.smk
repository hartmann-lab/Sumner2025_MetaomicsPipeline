
rule get_manifest:
    output:
        RESULTS + "AMP/qiime2/setup/manifest.txt"
    run:
        with open(output[0], "w+") as f:
            f.write("sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n")
            for sample in samples_AMP["sample"]:
                f.write(f"{samples.loc[sample, "sample"]}\t{samples.loc[sample, "AMP_read1"]}\t{samples.loc[sample, "AMP_read2"]}\n")


rule qiime_import:
    input:
        manifest = RESULTS + "AMP/qiime2/setup/manifest.txt",
        reads1 = samples_AMP['AMP_read1'], # Check reads actually there
        reads2 = samples_AMP['AMP_read2']
    output:
        RESULTS + "AMP/qiime2/artifacts/amp.qza"
    shell:
        """
        module load qiime2/2021.11
        qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path {input.manifest} --output-path {output} --input-format PairedEndFastqManifestPhred33V2
        """


rule dada2:
    input:
        amp = RESULTS + "AMP/qiime2/artifacts/amp.qza"
    output:
        rep_seqs = RESULTS + "AMP/qiime2/artifacts/dada2.qza",
        table = RESULTS + "AMP/qiime2/artifacts/table.qza",
        stats = RESULTS + "AMP/qiime2/artifacts/stats.qza"
    params:
        trim_left_f = 17,
        trim_left_r = 21,
        trunc_left_f = 270,
        trunc_left_r =270
    threads: 25
    shell:
        """
        module load qiime2/2021.11
        qiime dada2 denoise-paired --i-demultiplexed-seqs {input.amp} --p-trim-left-f {params.trim_left_f} --p-trim-left-r {params.trim_left_r} --p-trunc-len-f {params.trunc_left_f} --p-trunc-len-r {params.trunc_left_r} --o-representative-sequences {output.rep_seqs} --o-table {output.table} --o-denoising-stats {output.stats} --p-no-hashed-feature-ids --p-n-threads {threads} --verbose
        """


rule summarize:
    input:
        rep_seqs = RESULTS + "AMP/qiime2/artifacts/dada2.qza",
        table = RESULTS + "AMP/qiime2/artifacts/table.qza",
        stats = RESULTS + "AMP/qiime2/artifacts/stats.qza",
        taxa = RESULTS + "AMP/qiime2/artifacts/taxonomy.qza"
    output:
        rep_seqs =RESULTS + "AMP/qiime2/visuals/dada2.qzv",
        table = RESULTS + "AMP/qiime2/visuals/table.qzv",
        stats = RESULTS + "AMP/qiime2/visuals/stats.qzv",
        taxa = RESULTS + "AMP/qiime2/visuals/taxonomy.qzv"
    shell:
        """
        module load qiime2/2021.11
        qiime feature-table summarize --i-table {input.table} --o-visualization {output.table}
        qiime feature-table tabulate-seqs --i-data {input.rep_seqs} --o-visualization {output.rep_seqs}
        qiime metadata tabulate --m-input-file {input.stats} --o-visualization {output.stats}
        qiime metadata tabulate --m-input-file {input.taxa} --o-visualization {output.taxa}
        """


rule phylogeny:
    input:
        rep_seqs = RESULTS + "AMP/qiime2/artifacts/dada2.qza"
    output:
        rooted_tree = RESULTS + "AMP/qiime2/artifacts/rooted_tree.qza",
        unrooted_tree = RESULTS + "AMP/qiime2/artifacts/unrooted_tree.qza",
        alignmnet = RESULTS + "AMP/qiime2/artifacts/alignment.qza",
        masked_alignment = RESULTS + "AMP/qiime2/artifacts/masked_alignment.qza"
    threads: 24
    shell:
        """
        module load qiime2/2021.11
        qiime phylogeny align-to-tree-mafft-fasttree --i-sequences {input.rep_seqs} --o-alignment {output.alignmnet} --o-masked-alignment {output.masked_alignment} --o-tree {output.unrooted_tree} --o-rooted-tree {output.rooted_tree} --p-n-threads {threads} --verbose
        """


rule get_qiime2_database:
    output:
        RESULTS + "databases/qiime2_silva/silva-138-99-seqs.qza",
        RESULTS + "databases/qiime2_silva/silva-138-99-tax.qza"
    shell:
        """
        cd $(dirname {output[0]}) 
        wget https://data.qiime2.org/2021.11/common/silva-138-99-seqs.qza
        wget https://data.qiime2.org/2021.11/common/silva-138-99-tax.qza 
        """


rule classifer:
    input:
        seqs = RESULTS + "databases/qiime2_silva/silva-138-99-seqs.qza",
        tax = RESULTS + "databases/qiime2_silva/silva-138-99-tax.qza",
        rep_seqs = RESULTS + "AMP/qiime2/artifacts/dada2.qza"
    output:
        RESULTS + "AMP/qiime2/artifacts/taxonomy.qza"
    shell:
        """
        module load qiime2/2021.11
        qiime feature-classifier classify-consensus-vsearch --i-query {input.rep_seqs} --i-reference-reads {input.seqs} --i-reference-taxonomy {input.tax} --o-classification {output} --p-threads {threads} --verbose
        """


rule barplot:
    input:
        table = RESULTS + "AMP/qiime2/artifacts/table.qza",
        taxonomy = RESULTS + "AMP/qiime2/artifacts/taxonomy.qza",
    output:
        RESULTS + "AMP/qiime2/visuals/barplot.qzv"
    shell:
        """
        module load qiime2/2021.11
        qiime taxa barplot --i-table {input.table} --i-taxonomy {input.taxonomy} --o-visualization {output} --verbose
        """