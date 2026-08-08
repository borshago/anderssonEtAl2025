SAMPLES = ["PK8207_Pop1","PK8207_Pop2","PK8208_Pop1","PK8208_Pop2","PK8209_Pop1","PK8209_Pop2","PK8293_Pop1","PK8293_Pop2","PK2323_omito", "PK2323_ymito", "PK2323_ymitoPCs", "PK2326_ymito", "PK2326_omito", "PK2326_ymitoPCs"]

rule all:
    """
    Collect the main outputs of the workflow.
    """
    input:
        "RESULTS/rulegraph.png",
        #"DATA/GENOME/GRCm38.primary_assembly.genome.fa",
        #"DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf",
        #expand("RESULTS/QCOUT/fastqc_out/{samples}_R1_fastqc.html",samples=SAMPLES),
        #expand("RESULTS/QCOUT/fastqc_out/{samples}_R2_fastqc.html",samples=SAMPLES),
        #expand("RESULTS/QCOUT/fastqc_out/{samples}_R3_fastqc.html",samples=SAMPLES)
        #"DATA/GENOME/STAR/Genome",
        #expand("DATA/PROCESSED/umi_{samples}_R1.fastq.gz",samples=SAMPLES),
        #expand("DATA/PROCESSED/umi_{samples}_R3.fastq.gz",samples=SAMPLES)
        #expand("DATA/ALIGNED/{samples}/{samples}_Aligned.sortedByCoord.out.bam",samples=SAMPLES),
        #expand("DATA/ALIGNED/{samples}/{samples}_Aligned.sortedByCoord.out.bam.bai",samples=SAMPLES)
        #expand("DATA/ALIGNED/{samples}/{samples}_dedup.aligned.bam",samples=SAMPLES),
        #expand("RESULTS/QCOUT/fastqc_out_dedup/{samples}_dedup.aligned_fastqc.html",samples=SAMPLES)
        #expand("DATA/ALIGNED/dedup_{samples}_R1.bam",samples=SAMPLES),
        #expand("DATA/ALIGNED/dedup_{samples}_R3.bam",samples=SAMPLES),
        #expand("DATA/PROCESSED/dedup_{samples}_R1.fastq",samples=SAMPLES),
        #expand("DATA/PROCESSED/dedup_{samples}_R3.fastq",samples=SAMPLES)
        #"DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed",
        #"DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12",
        #expand("RESULTS/QCOUT/RSEQC/strandedness_{samples}.txt",samples=SAMPLES),
        #expand("RESULTS/QCOUT/RSEQC/rpkmSaturation_{samples}.saturation.pdf",samples=SAMPLES),
        #expand("RESULTS/QCOUT/RSEQC/junctionSaturation_{samples}.junctionSaturation_plot.pdf",samples=SAMPLES),
        #"DATA/GENOME/gencode.vM24.transcripts.idx",
        #expand("DATA/COUNTS/{samples}/abundance.tsv",samples=SAMPLES),
        #expand("RESULTS/QCOUT/RSEQC/readDistribution_{samples}.txt",samples=SAMPLES),
        #expand("DATA/ALIGNED/{samples}/out/aligned.log",samples=SAMPLES),
        #"RESULTS/OLD_VS_YOUNG//de_gsea_analysis.html",
        #"RESULTS/PANETH_CONTAM/de_gsea_panethContam.html"

rule downloadGenome:
    """
    Download GENCODE mouse genome and annotations.
    """
    output:
        "DATA/GENOME/GRCm38.primary_assembly.genome.fa",
        "DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf",
        "DATA/GENOME/gencode.vM24.transcripts.fa.gz"
    shell:
        """
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/GRCm38.primary_assembly.genome.fa.gz -P DATA/GENOME
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.primary_assembly.annotation.gtf.gz -P DATA/GENOME
        wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.transcripts.fa.gz -P DATA/GENOME
        gunzip DATA/GENOME/GRCm38.primary_assembly.genome.fa.gz
        gunzip DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf.gz
        """

rule annotationToBed:
    """
    Convert GTF annotation into BED format.
    """
    input:
        "DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf"
    output:
        "DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12"
    shell:
        """
        SOFTWARE/gtfToGenePred DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf DATA/GENOME/gencode.vM24.primary_assembly.annotation.genePred
        SOFTWARE/genePredToBed DATA/GENOME/gencode.vM24.primary_assembly.annotation.genePred DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12
        """

rule fastqc:
    """
    Run FastQC on a FASTQ file.
    """
    input:
        "DATA/READS/{sample}_R1.fastq.gz",
        "DATA/READS/{sample}_R2.fastq.gz",
        "DATA/READS/{sample}_R3.fastq.gz"
    output:
        "RESULTS/QCOUT/fastqc_out/{sample}_R1_fastqc.html",
        "RESULTS/QCOUT/fastqc_out/{sample}_R2_fastqc.html",
        "RESULTS/QCOUT/fastqc_out/{sample}_R3_fastqc.html"
    shell:
        """
        # Run fastQC and save the output to the QCOUT directory
        fastqc {input} -o RESULTS/QCOUT/fastqc_out
        """

rule createIndex:
    """
    Create STAR genome index.
    """
    input:
        "DATA/GENOME/GRCm38.primary_assembly.genome.fa",
        "DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf"
    output:
        "DATA/GENOME/STAR/Genome"
    shell:
        """
        STAR --runThreadN 20 --runMode genomeGenerate --genomeDir DATA/GENOME/STAR --genomeFastaFiles DATA/GENOME/GRCm38.primary_assembly.genome.fa --sjdbGTFfile DATA/GENOME/gencode.vM24.primary_assembly.annotation.gtf
        """

rule extractUmis:
    """
    Cut UMI sequence from Read 2 and add it to the name of Read 1 and 3
    """
    input:
        "DATA/READS/{sample}_R1.fastq.gz",
        "DATA/READS/{sample}_R2.fastq.gz",
        "DATA/READS/{sample}_R3.fastq.gz"
    output:
        "DATA/PROCESSED/umi_{sample}_R1.fastq.gz",           
        "DATA/PROCESSED/umi_{sample}_R3.fastq.gz"
    log:
        "LOGS/extractUmis_{sample}.log.txt"
    shell:
        """
        umi_tools extract -I {input[1]} --bc-pattern=NNNNNNNN --read2-in={input[0]} --stdout=/dev/null --read2-out={output[0]} >> {log}
        umi_tools extract -I {input[1]} --bc-pattern=NNNNNNNN --read2-in={input[2]} --stdout=/dev/null --read2-out={output[1]} >> {log}
        """

rule map:
    """
    Map UMI-extracted reads.
    """
    input:
        "DATA/PROCESSED/umi_{sample}_R1.fastq.gz",
        "DATA/PROCESSED/umi_{sample}_R3.fastq.gz"
    output:
        "DATA/ALIGNED/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    threads: 2
    shell:
        """
        STAR --readFilesIn {input[0]} {input[1]} --runThreadN {threads} --genomeDir DATA/GENOME/STAR --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix DATA/ALIGNED/{wildcards.sample}/{wildcards.sample}_ --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0.5 --outFilterMatchNminOverLread 0.5
        """

rule indexBam:
    """
    Create bam indices.
    """
    input:
        "DATA/ALIGNED/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "DATA/ALIGNED/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        """
        samtools index {input} {output}
        """

rule deduplicate:
    """
    Deduplicate the mapped reads using UMIs.
    """
    input:
        "DATA/ALIGNED/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        "DATA/ALIGNED/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam",
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam.bai"
    log:
        "LOGS/deduplicate_{sample}.log.txt"
    shell:
        """
        umi_tools dedup -I {input[0]} --paired --output-stats=DATA/ALIGNED/{wildcards.sample} -S {output[0]} > {log}
        samtools index {output} # test this!
        """

rule fastqcDedup:
    """
    Run FastQC on a deduplicated FASTQ file.
    """
    input:
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam"
    output:
        "RESULTS/QCOUT/fastqc_out_dedup/{sample}_dedup.aligned_fastqc.html"
    shell:
        """
        # Run fastQC and save the output to the QCOUT directory
        fastqc {input} -o RESULTS/QCOUT/fastqc_out_dedup
        """

rule splitBam:
    """
    Split aligned and deduplicated reads for quantification with kallisto.
    """
    input:
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam"
    output:
        "DATA/ALIGNED/dedup_{sample}_R1.bam",
        "DATA/ALIGNED/dedup_{sample}_R3.bam"
    shell:
        """
        samtools view -b -h -f 0x1 -f 0x40 DATA/ALIGNED/{wildcards.sample}/{wildcards.sample}_dedup.aligned.bam > DATA/ALIGNED/dedup_{wildcards.sample}_R1.bam
        samtools view -b -h -f 0x1 -f 0x80 DATA/ALIGNED/{wildcards.sample}/{wildcards.sample}_dedup.aligned.bam > DATA/ALIGNED/dedup_{wildcards.sample}_R3.bam
        """

rule bamtoFastq:
    """
    Convert split BAM files to FASTQ.
    """
    input:
        "DATA/ALIGNED/dedup_{sample}_R1.bam",
        "DATA/ALIGNED/dedup_{sample}_R3.bam"
    output:
        "DATA/PROCESSED/dedup_{sample}_R1.fastq",
        "DATA/PROCESSED/dedup_{sample}_R3.fastq"
    shell:
        """
        samtools bam2fq DATA/ALIGNED/dedup_{wildcards.sample}_R1.bam > DATA/PROCESSED/dedup_{wildcards.sample}_R1.fastq
        samtools bam2fq DATA/ALIGNED/dedup_{wildcards.sample}_R3.bam > DATA/PROCESSED/dedup_{wildcards.sample}_R3.fastq
        """

rule rseqc:
    """
    Run post-alignment QC.
    """
    input:
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam",
        "DATA/ALIGNED/{sample}/{sample}_dedup.aligned.bam.bai"
    output:
        "RESULTS/QCOUT/RSEQC/strandedness_{sample}.txt",
        "RESULTS/QCOUT/RSEQC/rpkmSaturation_{sample}.saturation.pdf",
        "RESULTS/QCOUT/RSEQC/junctionSaturation_{sample}.junctionSaturation_plot.pdf",
        "RESULTS/QCOUT/RSEQC/readDistribution_{sample}.txt"
    shell:
        """
        infer_experiment.py -r DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12 -i {input[0]} > {output[0]}
        geneBody_coverage.py -r DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12 -i {input[0]} -o RESULTS/QCOUT/RSEQC/coverage_{wildcards.sample}
        RPKM_saturation.py -r DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12 -d '1++,1--,2+-,2-+' -q 255 -i {input[0]} -o RESULTS/QCOUT/RSEQC/rpkmSaturation_{wildcards.sample}
        junction_saturation.py -i {input[0]} -r DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12 -o RESULTS/QCOUT/RSEQC/junctionSaturation_{wildcards.sample}
        read_distribution.py  -i {input[0]} -r DATA/GENOME/gencode.vM24.primary_assembly.annotation.bed12 > RESULTS/QCOUT/RSEQC/readDistribution_{wildcards.sample}.txt
        """

rule kallistoIndex:
    """
    Create kallisto index.
    """
    input:
        "DATA/GENOME/gencode.vM24.transcripts.fa.gz"
    output:
        "DATA/GENOME/gencode.vM24.transcripts.idx"
    shell:
        """
        kallisto index -i DATA/GENOME/gencode.vM24.transcripts.idx DATA/GENOME/gencode.vM24.transcripts.fa.gz
        """

rule kallistoQuant:
    """
    Quantify gene expression with kallisto.
    """
    input:
        "DATA/PROCESSED/dedup_{sample}_R1.fastq",
        "DATA/PROCESSED/dedup_{sample}_R3.fastq"
    output:
        "DATA/COUNTS/{sample}/abundance.tsv"
    shell:
        """
        kallisto quant -t 40 --plaintext --fr-stranded -i DATA/GENOME/gencode.vM24.transcripts.idx -o DATA/COUNTS/{wildcards.sample} {input[0]} {input[1]}
        """

rule deGseaAnalysis:
    """
    Run R scripts for differential expression and gene set enrichment analyses.
    """
    input:
        #expand("DATA/COUNTS/{samples}/abundance.tsv", samples=SAMPLES),
        de = "CODE/de_gsea_analysis.R",
        contam = "CODE/de_gsea_panethContam.R"
    output:
        "RESULTS/OLD_VS_YOUNG/de_gsea_analysis.html",
        "RESULTS/PANETH_CONTAM/de_gsea_panethContam.html"
    shell:
        """
        mkdir RESULTS/OLD_VS_YOUNG RESULTS/OLD_VS_YOUNG/DE RESULTS/OLD_VS_YOUNG/GSEA RESULTS/PANETH_CONTAM RESULTS/PANETH_CONTAM/DE RESULTS/PANETH_CONTAM/GSEA 
        Rscript -e 'rmarkdown::render({input.de}, output_file={output[0]})'
        Rscript -e 'rmarkdown::render({input.contam}, output_file={output[1]})'
        """

rule generate_rulegraph:
    """
    Generate a rulegraph for the workflow.
    """
    output:
        "RESULTS/rulegraph.png"
    shell:
        """
        snakemake --snakefile Snakefile --rulegraph | dot -Tpng > {output}
        """
