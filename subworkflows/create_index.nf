include { INDEX_REF; FAI_REF; DICT_REF } from '../modules/modules.nf'

workflow CREATE_INDEX {
    take: 
        reference

    main:
        INDEX_REF(reference)
        FAI_REF(reference)
        DICT_REF(reference)

    emit:
        bwa_index = INDEX_REF.out
        fai_index = FAI_REF.out
        seq_dict = DICT_REF.out
}