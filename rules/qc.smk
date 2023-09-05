
if config['options']['paired']:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}_{pair}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}_{pair}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}_{pair}_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_raw')
        conda:
            '../envs/fastqc.yaml'
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'

    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}_{pair}_trimmed.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}_{pair}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}_{pair}_trimmed_fastqc.zip')
        params:
            output_dir = path.join(qc_dir, 'fastqc_trim')
        conda:
            '../envs/fastqc.yaml'
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


else:

    rule fastqc_raw:
        input:
            path.join(raw_dir, '{sample}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_raw', '{sample}_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_raw', '{sample}_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_raw')
        conda:
            '../envs/fastqc.yaml'
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


    rule fastqc_trim:
        input:
            path.join(raw_dir, '{sample}.fastq.gz')
        output:
            html = path.join(qc_dir, 'fastqc_trim', '{sample}_trimmed_fastqc.html'),
            zip = path.join(qc_dir, 'fastqc_trim', '{sample}_trimmed_fastqc.zip')
        params:
            output_dir=path.join(qc_dir, 'fastqc_trim')
        conda:
            '../envs/fastqc.yaml'
        shell:
            'mkdir -p {params.output_dir} && '
            'fastqc --quiet --outdir {params.output_dir} {input}'


rule samtools_stats:
    input:
        path.join(bam_dir, '{sample}','Aligned.sortedByCoord.out.bam')
    output:
        path.join(stats_dir, '{sample}.txt')
    conda:
        '../envs/samtools.yaml'
    shell:
        'samtools stats {input} > {output}'


def multiqc_prerequisite(wildcards):

# fastqc for raw files
    fastqc_temp_raw = expand(path.join(qc_dir, 'fastqc_raw', '{sample}{{pair}}_fastqc.zip'), sample = all_samples)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_raw = expand(fastqc_temp_raw, pair=pairs)

# fastqc for trimmed file
    fastqc_temp_trim = expand(path.join(qc_dir, 'fastqc_trim', '{sample}{{pair}}_trimmed_fastqc.zip'), sample = all_samples)
    pairs = ["_R1", "_R2"] if config['options']['paired'] else [""]
    fastqc_trim = expand(fastqc_temp_trim, pair=pairs)

# samtools stats
    samtools_stats = expand(path.join(stats_dir, '{sample}.txt'), sample = all_samples)

# quantification results
    quant_out = path.join(featureCounts_dir, 'merged.gene.txt')
    
    return fastqc_raw + fastqc_trim + samtools_stats + [quant_out]

#rule multiqc:
#    input:
#        multiqc_prerequisite
#    output:
#        path.join(qc_dir, 'multiqc_report.html')
#    params:
#        qc_dir = qc_dir,
#        out_dir = qc_dir
#    conda:
#        '../envs/multiqc.yaml'
#    log:
#        path.join(log_dir, 'multiqc.log')
#    shell:
#        'multiqc --outdir {params.out_dir} {params.qc_dir} 2> {log}'

rule qc_merge:
    input:
        multiqc_prerequisite
    output: 
        touch(path.join(qc_dir, 'qc_files.done'))


