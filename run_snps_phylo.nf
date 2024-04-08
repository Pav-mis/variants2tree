include { CREATE_INDEX} from "./subworkflows/create_index.nf"
include { UBAM_QC_AND_MAPPING} from "./subworkflows/ubam_qc_and_mapping.nf"
include { VCF_GENOTYPING_AND_FILTERING} from "./subworkflows/VCF_genotyping_and_filtering.nf"
include { BUILD_TREE} from "./subworkflows/build_tree.nf"
include { GENERATE_UBAM; VARIANT_CALLING} from "./modules/modules.nf"

workflow RUN_SNPS_PHYLO {
    reference = file(params.ref) //define input channel
    fastq = Channel.fromFilePairs("${params.fastq}/*_{R1,R2}.fastq") //define input channel
    CREATE_INDEX(reference) //subworkflow
    GENERATE_UBAM(fastq) //module
    UBAM_QC_AND_MAPPING(GENERATE_UBAM.out, reference, CREATE_INDEX.out.bwa_index, CREATE_INDEX.out.seq_dict) //subworkflow
    VARIANT_CALLING(UBAM_QC_AND_MAPPING.out, reference, CREATE_INDEX.out.fai_index, CREATE_INDEX.out.seq_dict) //module
    VCF_GENOTYPING_AND_FILTERING(VARIANT_CALLING.out.collect(), reference, CREATE_INDEX.out.fai_index, CREATE_INDEX.out.seq_dict) //subworkflow
    BUILD_TREE(CF_GENOTYPING_AND_FILTERING.out) //subworkflow
    
}