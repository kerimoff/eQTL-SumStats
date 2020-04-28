SNP_DSET = 'variant'
MANTISSA_DSET = 'mantissa'
EXP_DSET = 'exponent'
PVAL_DSET = 'pvalue'
STUDY_DSET = 'study_id'
TISSUE_DSET = 'tissue' # ontology term
CHR_DSET = 'chromosome'
BP_DSET = 'position'
BETA_DSET = 'beta'
SE_DSET = 'se'
EFFECT_DSET = 'alt'
OTHER_DSET = 'ref'
FREQ_DSET = 'maf'
RSID_DSET = 'rsid'
MUTATION_DSET = 'type'
PHEN_DSET = 'molecular_trait_id'
AC_DSET = 'ac'
AN_DSET = 'an'
R2_DSET = 'r2'
GENE_DSET = 'gene_id'
MTO_DSET = 'molecular_trait_object_id'
EXPR_DSET = 'median_tpm'
HM_OR_DSET = 'hm_odds_ratio'
HM_RANGE_U_DSET = 'hm_ci_upper'
HM_RANGE_L_DSET = 'hm_ci_lower'
HM_BETA_DSET = 'hm_beta'
HM_EFFECT_DSET = 'hm_effect_allele'
HM_OTHER_DSET = 'hm_other_allele'
HM_FREQ_DSET = 'hm_effect_allele_frequency'
HM_VAR_ID = 'hm_variant_id'
HM_CODE = 'hm_code'
QTL_GROUP_DSET = 'qtl_group'
CONDITION_DSET = 'condition'
CONDITION_LABEL_DSET = 'condition_label'
TISSUE_LABEL_DSET = 'tissue_label'

#qtl_group, condition, condition_label, cell_type, ontology_term, ontology_label


DSET_TYPES = {SNP_DSET: str, RSID_DSET: str, MUTATION_DSET: str, AC_DSET: float, AN_DSET: float, PVAL_DSET: float, MANTISSA_DSET: float, EXP_DSET: "int64", STUDY_DSET: str,
              CHR_DSET: str, BP_DSET: "int64", R2_DSET: float, BETA_DSET: float, SE_DSET: float,
              EFFECT_DSET: str, OTHER_DSET: str, FREQ_DSET: float, EXPR_DSET: float, TISSUE_DSET: str,
              QTL_GROUP_DSET: str, CONDITION_DSET: str, CONDITION_LABEL_DSET: str, TISSUE_LABEL_DSET: str}
              

REFERENCE_DSET = SNP_DSET
HARMONISATION_PREFIX = 'hm_'
GWAS_CATALOG_STUDY_PREFIX = 'GCST'

TO_DISPLAY_DEFAULT = {SNP_DSET, PVAL_DSET, STUDY_DSET, CHR_DSET, BP_DSET, EFFECT_DSET, OTHER_DSET, BETA_DSET, RSID_DSET, MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET, R2_DSET, EXPR_DSET, QTL_GROUP_DSET, CONDITION_DSET, CONDITION_LABEL_DSET, TISSUE_LABEL_DSET}

TO_DISPLAY_RAW = {SNP_DSET, PVAL_DSET, STUDY_DSET, CHR_DSET, BP_DSET, BETA_DSET,
                  EFFECT_DSET, OTHER_DSET}


TO_LOAD_DSET_HEADERS_DEFAULT = {PHEN_DSET, SNP_DSET, PVAL_DSET, CHR_DSET, BP_DSET, EFFECT_DSET, OTHER_DSET, BETA_DSET,  MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET, R2_DSET, EXPR_DSET, GENE_DSET, MTO_DSET, RSID_DSET, SE_DSET}
#TO_STORE_DSETS_DEFAULT = {SNP_DSET, MANTISSA_DSET, EXP_DSET, STUDY_DSET, CHR_DSET, BP_DSET, EFFECT_DSET, OTHER_DSET, BETA_DSET, RSID_DSET, MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET,  R2_DSET, EXPR_DSET}
#TO_QUERY_DSETS_DEFAULT = {SNP_DSET, MANTISSA_DSET, EXP_DSET, STUDY_DSET, CHR_DSET, BP_DSET, BETA_DSET, RSID_DSET, MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET, R2_DSET, MEAN_EXPR_DSET,
#                  EFFECT_DSET, OTHER_DSET}
# temp change tp pvalue instead of mantissa exp.
TO_STORE_DSETS_DEFAULT = {SNP_DSET, PVAL_DSET, STUDY_DSET, CHR_DSET, BP_DSET, EFFECT_DSET, OTHER_DSET, BETA_DSET, RSID_DSET, MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET, SE_DSET, R2_DSET, EXPR_DSET}
TO_QUERY_DSETS_DEFAULT = {SNP_DSET, PVAL_DSET, STUDY_DSET, CHR_DSET, BP_DSET, BETA_DSET, RSID_DSET, MUTATION_DSET, AC_DSET, AN_DSET, FREQ_DSET, R2_DSET, EXPR_DSET,
                          EFFECT_DSET, OTHER_DSET, TISSUE_DSET}
TO_INDEX = [PHEN_DSET, BP_DSET, PVAL_DSET, SNP_DSET, GENE_DSET]
TRAIT_FILE_INDEX = ['phenotype_id', 'gene_id']
CHROMOSOMES = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
