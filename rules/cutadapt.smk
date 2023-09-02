################################################################################ 
# trimming by cutadapt
################################################################################ 
# If raw file is bam foramt, it is transformed into fastq format
# Before that, the reads in the bam file should be sorted according to its read name
# This is for cutadapt which requires R1-R2 read matches for paired data

if config['options']['paired']:

    rule cutadapt:
        input:
            R1 = path.join(raw_dir, '{sample}_R1.fastq.gz'),
            R2 = path.join(raw_dir, '{sample}_R2.fastq.gz')
        output:
            R1 = temp(path.join(raw_dir, '{sample}_R1_trimmed.fastq.gz')),
            R2 = temp(path.join(raw_dir, '{sample}_R2_trimmed.fastq.gz'))
        params:
            options = format_options(config['cutadapt']['options'])
        conda:
            '../envs/cutadapt.yaml'
        threads:
            config['cutadapt']['threads']
        log:
            path.join(log_dir, 'cutadapt', '{sample}.log')
        shell:
            'cutadapt {params.options} --cores={threads} -o {output.R1} -p {output.R2} '
            '{input.R1} {input.R2} 2> {log}' 

else:
   
    rule cutadapt:
        input:
            path.join(raw_dir, '{sample}.fastq.gz')
        output:
            temp(path.join(raw_dir, '{sample}_trimmed.fastq.gz'))
        params:
            options = format_options(config['cutadapt']['options'])
        conda:
            '../envs/cutadapt.yaml'
        threads:
            config['cutadapt']['threads']
        log:
            path.join(log_dir, 'cutadapt', '{sample}.log'),
        shell:
            'cutadapt {params.options} --cores={threads} -o {output} {input} 2> {log}'



