
def get_reads(wildcards):
    return samples.loc[wildcards.sample, ["sample", "AMP_read1", "AMP_read2", "MGX_read1", "MGX_read2", "MTX_read1", "MTX_read2"]].dropna()

def get_reads_mgx(wildcards):
    return samples_MGX.loc[wildcards.sample, ["sample", "AMP_read1", "AMP_read2", "MGX_read1", "MGX_read2", "MTX_read1", "MTX_read2"]].dropna()

def get_amp_read1(wildcards):
    return get_reads(wildcards)["AMP_read1"]

def get_amp_read2(wildcards):
    return get_reads(wildcards)["AMP_read2"]

def get_mgx_read1(wildcards):
    return get_reads_mgx(wildcards)["MGX_read1"]

def get_mgx_read2(wildcards):
    return get_reads_mgx(wildcards)["MGX_read2"]

def get_mtx_read1(wildcards):
    return get_reads(wildcards)["MTX_read1"]
