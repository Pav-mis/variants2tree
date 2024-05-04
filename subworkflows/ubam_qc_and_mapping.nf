include { SAM_TO_FASTQ; MARK_ILLUMINA_ADAPTERS; ALIGN_TO_REF; SORT_AND_INDEX_BAM; MARK_DUPLICATES; MERGE_BAM_WITH_UBAM} from '../modules/modules.nf'

workflow UBAM_QC_AND_MAPPING {
    take:
        reads
        ubams
        reference
        bwa_index
        seq_dict

    main:
        MARK_ILLUMINA_ADAPTERS(reads, ubams)
        ALIGN_TO_REF(MARK_ILLUMINA_ADAPTERS.out, reference, bwa_index)
        SORT_AND_INDEX_BAM(ALIGN_TO_REF.out)
        MARK_DUPLICATES(SORT_AND_INDEX_BAM.out)
        MERGE_BAM_WITH_UBAM(MARK_DUPLICATES.out, reference, seq_dict)
    
    emit:
        ubam = MERGE_BAM_WITH_UBAM.out[0]
        bam = MERGE_BAM_WITH_UBAM.out[1]
        bam_index = MERGE_BAM_WITH_UBAM.out[2]
        
    
}


    