
rule move_tables_to_final_mgx:
    input:
        qc_counts = 
        csome_counts = 
        mpa = 
        humann = 
    output: 
        qc_counts = 
        csome_counts = 
        mpa = 
        humann = 
    params:
        out_dir = 
    shell:
        """
        cp {input} {params.out_dir}
        """


rule move_directories:
    input:
        mgx = <MGX final>
        mtx = <MTX final>
    output:
        mgx = <MGX final>
        mtx = <MTX final>
    params:
        mgx_dir = RESULTS + "output/mgx/"
        mtx_dir =  RESULTS + "output/mtx/"
    shell:
        """
        cp -R {input.mgx} {params.mgx_dir}
        cp -R {input.mtx} {params.mtx_dir}
        """