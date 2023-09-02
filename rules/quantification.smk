# FeatureCounts for raw read counts
rule featureCounts_gene:
    input:
        bam = path.join(bam_dir, '{sample}', 'Aligned.sortedByCoord.out.bam'),
        bai = path.join(bam_dir, '{sample}', 'Aligned.sortedByCoord.out.bam.bai'),
        star_unload = path.join(bam_dir, 'star_unload.done')
    output:
        path.join(featureCounts_dir, '{sample}_counts.gene.txt')
    params:
        options = format_options(config['featureCounts']['options']),
        gtf = config['featureCounts']['gtf']
    conda:
            '../envs/featurecounts.yaml'
    threads:
        config['featureCounts']['threads']
    log:
        path.join(log_dir, 'featureCounts', '{sample}_gene.txt')
    shell:
        'featureCounts {params.options} -a {params.gtf} '
        '-t exon -g gene_id -o {output} -T {threads} {input.bam} 2> {log}'

rule merge_gene:
    input:
        expand(path.join(featureCounts_dir, '{sample}_counts.gene.txt'),
            sample = all_samples)
    output:
        path.join(featureCounts_dir, 'merged.gene.txt')
    run:
        # Merge count files.
        frames = (pd.read_csv(fp, sep="\t", skiprows=1, index_col=list(range(6)))
                    for fp in input)
        merged = pd.concat(frames, axis=1)
        # Extract sample names
        merged = merged.rename(
                    columns=lambda c: path.splitext(path.basename(c))[0])
        
        merged.to_csv(output[0], sep="\t", index=True)

