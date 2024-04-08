include { COMBINE_AND_GENOTYPE_VCF; FILTER_SNPS_AND_INDELS; QUALITY_FILTER_VARIANTS; FINAL_FILTER_VARIANTS} from '../modules/modules.nf'

workflow VCF_GENOTYPING_AND_FILTERING {
    take:
        gvcfs
        reference
        fai_index
        seq_dict

    main:
        COMBINE_AND_GENOTYPE_VCF(gvcfs, reference, fai_index, seq_dict)
        FILTER_SNPS_AND_INDELS(COMBINE_AND_GENOTYPE_VCF.out, reference, fai_index, seq_dict)
        QUALITY_FILTER_VARIANTS(FILTER_SNPS_AND_INDELS.out, reference, fai_index, seq_dict)
        FINAL_FILTER_VARIANTS(QUALITY_FILTER_VARIANTS.out)
    
    emit:
        FINAL_FILTER_VARIANTS.out
}