include { VCF_TO_PHYLIP; GENERATE_TREE} from '../modules/modules.nf'

workflow BUILD_TREE {
    take:
        vcf

    main:
        VCF_TO_PHYLIP(vcf)
        GENERATE_TREE(VCF_TO_PHYLIP.out)
    
    emit:
        tree = GENERATE_TREE.out
}