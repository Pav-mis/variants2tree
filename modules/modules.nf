process INDEX_REF {
    input:
    path reference 

    output:
    path "${reference}*"

    """
    bwa index -a bwtsw ${reference}
    """
}

process FAI_REF {
    input:
    path reference 

    output:
    path "${reference}.fai"

    """
    samtools faidx ${reference}
    """
}

process DICT_REF {
    input:
    path reference 

    output:
    path "*.dict"

    """
    gatk CreateSequenceDictionary -R ${reference}
    """
}

process GENERATE_UBAM {
    errorStrategy 'ignore'
    
    input:
    tuple val(sampleId), file(reads)

    output:
    tuple val(sampleId), file(reads)
    path "${sampleId}.ubam"

    """
    mkdir tmp
    platform="Illumina"
    seqprovider="AGRF"
    date="2023-01-40T00:00:00-0400"
    gatk FastqToSam \
        -F1 ${reads[0]} \
        -F2 ${reads[1]} \
        -O ${sampleId}".ubam" \
        --READ_GROUP_NAME ${sampleId} \
        --SAMPLE_NAME ${sampleId}\
        --LIBRARY_NAME ${sampleId} \
        --PLATFORM \$platform \
        --SEQUENCING_CENTER \$seqprovider \
        --RUN_DATE \$date \
        --TMP_DIR tmp
    """

}

process SAM_TO_FASTQ {
    input:
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${sampleId}.fastq"
    

    """
    gatk SamToFastq \
	    -I ${ubam} \
	    -FASTQ ${sampleId}.fastq \
	    -CLIPPING_ATTRIBUTE XT \
	    -CLIPPING_ACTION 2 \
	    -INTERLEAVE true \
	    -NON_PF true
    """
}

process MARK_ILLUMINA_ADAPTERS {
    input:
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${sampleId}.marked.fastq"
    tuple val(sampleId), file(reads)
    path ubam

    """
    mkdir tmp
    gatk MarkIlluminaAdapters \
        -I ${ubam} \
        -O ${sampleId}.marked.ubam \
        -M ${sampleId}.markilluminaadapters_metrics.txt \
        -TMP_DIR tmp
    gatk SamToFastq \
		-I ${sampleId}.marked.ubam \
		-FASTQ ${sampleId}.marked.fastq \
		-CLIPPING_ATTRIBUTE XT \
		-CLIPPING_ACTION 2 \
		-INTERLEAVE true \
		-NON_PF true 
    """
}

process ALIGN_TO_REF {
    input:
    path fastq
    tuple val(sampleId), file(reads)
    path ubam
    path reference
    path reference_index

    output:
    path "${fastq.simpleName}.sam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    bwa mem \
        -M \
        -t ${task.cpus} \
        -p ${reference} \
        ${fastq} > ${fastq.simpleName}.sam
    """
}

process SORT_AND_INDEX_BAM {
    input:
    path sam
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${sam.baseName}_v_${params.refname}.bam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    samtools sort ${sam} \
		-O BAM | \
		tee ${sam.baseName}_v_${params.refname}.bam | \
		samtools index - ${sam.baseName}_v_${params.refname}.bam.bai
    """
}

process MARK_DUPLICATES {
    input:
    path bam
    tuple val(sampleId), file(reads)
    path ubam

    output:
    path "${bam.baseName}.marked_duplicates.bam"
    tuple val(sampleId), file(reads)
    path ubam

    """
    gatk MarkDuplicates \
	    I=${bam} \
	    O=${bam.baseName}.marked_duplicates.bam \
	    M=${bam.baseName}.marked_duplicates.txt
    """
}

process MERGE_BAM_WITH_UBAM {
    input:
    path bam
    tuple val(sampleId), file(reads)
    path ubam
    path reference
    path seq_dict

    output:
    path ubam
    path "${bam.simpleName}.final.bam"
    path "${bam.simpleName}.final.bai"

    """
    echo ${ubam}
    mkdir tmp
    gatk MergeBamAlignment \
	    --ALIGNED ${bam} \
        --UNMAPPED ${ubam} \
	    -O ${bam.simpleName}.final.bam \
	    -R ${reference}
    gatk BuildBamIndex \
	    -I ${bam.simpleName}.final.bam \
	    -O ${bam.simpleName}.final.bai
    """
}

process VARIANT_CALLING {
    input:
    path ubam
    path bam
    path bamindex
    path reference
    path fai_index
    path seq_dict

    output:
    path "${bam.simpleName}.g.vcf" 

    """
    mkdir tmp
    gatk HaplotypeCaller \
        -I ${bam} \
        -O ${bam.simpleName}.g.vcf \
        -R ${reference} \
        -ERC GVCF \
        --minimum-mapping-quality 20 \
        --min-base-quality-score 20 \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        --tmp-dir tmp \
        --native-pair-hmm-threads ${task.cpus * 2}
    """ 
}

process COMBINE_AND_GENOTYPE_VCF {
    input:
    path gvcfs
    path reference
    path fai_index
    path seq_dict

    output:
    path "${params.refname}.combined_panel.vcf"
    path "${params.refname}.combined_panel.vcf.idx"

    """
    ls *.g.vcf > list.gvcfs.list
    gatk CombineGVCFs \
	    -R ${reference} \
	    --variant list.gvcfs.list \
	    -O ${params.refname}.combined_panel.g.vcf
    gatk GenotypeGVCFs \
	    -R ${reference} \
	    --variant ${params.refname}.combined_panel.g.vcf \
	    -O ${params.refname}.combined_panel.vcf
    """

}


process FILTER_SNPS_AND_INDELS {
    input:
    path vcf
    path vcf_index
    path reference
    path fai_index
    path seq_dict

    output:
    path "${params.refname}.combined_panel.snps.vcf"
    path "${params.refname}.combined_panel.snps.vcf.idx"
    path "${params.refname}.combined_panel.indels.vcf"
    path "${params.refname}.combined_panel.indels.vcf.idx"

    """
    gatk SelectVariants \
        -R ${reference} \
        -V ${vcf} \
        --select-type-to-include SNP \
        --output ${params.refname}.combined_panel.snps.vcf
    gatk SelectVariants \
        -R ${reference} \
        -V ${vcf} \
        --select-type-to-include INDEL \
        --output ${params.refname}.combined_panel.indels.vcf
    """

}

process QUALITY_FILTER_VARIANTS {
    input:
    path snp_vcf
    path snp_vcf_index
    path indel_vcf
    path indel_vcf_index
    path reference
    path fai_index
    path seq_dict

    output:
    path "${params.refname}.combined_panel.filtered.snps.vcf"
    path "${params.refname}.combined_panel.filtered.snps.vcf.idx"
    path "${params.refname}.combined_panel.filtered.indels.vcf"
    path "${params.refname}.combined_panel.filtered.indels.vcf.idx"

    """
    gatk VariantFiltration \
        -R ${reference} \
        -V ${snp_vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O ${params.refname}.combined_panel.filtered.snps.vcf
	
    gatk VariantFiltration \
        -R ${reference} \
        -V ${indel_vcf} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O ${params.refname}.combined_panel.filtered.indels.vcf
    """
}

process FINAL_FILTER_VARIANTS {
    input:
    path snps
    path snps_idx
    path indels
    path indels_idx

    output:
    path "${params.refname}.combined_panel.filtered.biallelic.1random5kb.ld90.snps.vcf"

    """
    bcftools view \
	--max-alleles 2 \
	--exclude-types indels \
	-e 'F_MISSING<=0.05' \
	${snps} -o ${params.refname}.combined_panel.filtered.biallelic.snps.bcf
	
    bcftools +prune \
        -m 0.9 \
        -w 5000bp \
        -n1 \
        -N rand \
        ${params.refname}.combined_panel.filtered.biallelic.snps.bcf \
        -Ov \
        -o ${params.refname}.combined_panel.filtered.biallelic.1random5kb.ld90.snps.vcf
    """
}

process VCF_TO_PHYLIP {
    input:
    path vcf

    output:
    path "${params.refname}.*.fasta"

    """
    python3 $projectDir/bin/vcf2phylip.py \
	    -i ${vcf} \
        -m 1 \
	    --fasta \
	    --nexus
    """
}

process GENERATE_TREE {
    input:
    path phylip

    output:
    path "*.treefile"

    """
    iqtree \
        -s ${phylip} \
        -st DNA \
        -m TEST \
        -bb 1000 \
        -alrt 1000
    """
}

