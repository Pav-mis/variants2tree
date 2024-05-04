include { INDEX_REF; FAI_REF; DICT_REF; FILTER_SNPS_AND_INDELS; QUALITY_FILTER_VARIANTS; FINAL_FILTER_VARIANTS; VCF_TO_PHYLIP; GENERATE_TREE} from "./modules/modules.nf"


workflow {
    reference = file(params.ref) //define input channel
    bwa_index = INDEX_REF(reference)
    fai_index = FAI_REF(reference)
    seq_dict = DICT_REF(reference)
    vcfs = file("/scratch/y95/pmisiun/example_results/A.fabae_59-13.fasta.combined_panel.vcf")
    vcf_idx = file("/scratch/y95/pmisiun/example_results/A.fabae_59-13.fasta.combined_panel.vcf.idx")
    snps_and_indels = FILTER_SNPS_AND_INDELS(vcfs, vcf_idx, reference, fai_index, seq_dict)
    quality_filtered_variants = QUALITY_FILTER_VARIANTS(snps_and_indels, reference, fai_index, seq_dict)
    final_variants = FINAL_FILTER_VARIANTS(quality_filtered_variants)
    phylip = VCF_TO_PHYLIP(final_variants)
    tree = GENERATE_TREE(phylip)
}