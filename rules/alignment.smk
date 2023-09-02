################################################################################ 
# STAR alignment
################################################################################ 

# input fastq files
def align_inputs(wildcards):
    fastq_path = path.join(raw_dir, "{sample}{{pair}}_trimmed.fastq.gz".format(
        sample=wildcards.sample))
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    return expand(fastq_path, pair = pairs)


# Using shared memory with genomeLoad
rule star_preload:
    params:
        index = config['star_align']['index']
    output:
        temp(touch(path.join(bam_dir, 'star_preload.done')))
    conda:
        '../envs/star.yaml'
    log:
        path.join(log_dir, 'star_align', 'genome_preload.log')
    shell:
        'STAR --genomeLoad LoadAndExit --genomeDir {params.index} 2> {log}'


rule star_align:
    input:
        fq = align_inputs,
        tmp = path.join(bam_dir, 'star_preload.done')
    output:
        protected(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'))
    params:
        index = config['star_align']['index'],
        options = format_options(config['star_align']['options']),
        out_prefix = path.join(bam_dir, '{sample}/'),  # '/' should be at the end of folder name: '{sample}/'
        qc_prefix = path.join(qc_dir,'star_align','{sample}'),
        qc_dir = path.join(qc_dir,'star_align')
    conda:
        '../envs/star.yaml'
    threads:
        config['star_align']['threads']
    log:
        path.join(log_dir, 'star_align', '{sample}.log')
    shell:
        'STAR {params.options} --genomeDir {params.index} '
        '--outFileNamePrefix {params.out_prefix} --runThreadN {threads} '
        '--readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate '
        '--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 '
        '--readFilesIn {input.fq} 2> {log} '
        '&& mkdir -p {params.qc_dir} '
        '&& mv {params.out_prefix}Log.final.out {params.qc_prefix}.Log.final.out'

# Genome unloading
rule star_unload:
    input:
        expand(path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam'),
            sample=all_samples)
    output:
        touch(path.join(bam_dir, 'star_unload.done'))
    params:
        index=config['star_align']['index']
    conda:
        '../envs/star.yaml'
    log:
        path.join(log_dir, 'star_align', 'genome_remove.log')
    shell:
        'STAR --genomeLoad Remove --genomeDir {params.index} > {log}'


################################################################################ 
# indexing the bam files
################################################################################ 

rule sambamba_index:
    input:
        path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam')
    output:
        path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam.bai')
    conda:
        '../envs/sambamba.yaml'
    threads:
        config['sambamba_index']['threads']
    shell:
        'sambamba index -t {threads} {input} {output}'


