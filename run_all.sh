#################################################
# Run all steps in LDSC regression for cell-type 
# specific heritability partitioning
##################################################

# Munge summary statistics
bash scripts/gwas_sumstat.sh

# Extract the necessary columns from the ATAC-seq data
bash scripts/prep_atac_seq.sh

# Lift over ATAC-seq peaks from hg38 to hg19
bash scripts/liftover_singlecell.sh

# Annotate SNPs by which tissues they belong to
Rscript scripts/tissue_specific_snp_annotation_single_cell.R

# Generate LD scores for each single cell type
bash scripts/ldscore_singlecell.sh

# Perform LD score regression to partition heritability
bash scripts/cell_type_spe_regression.sh

# Plot final enrichments
Rscript scripts/plot_cellspecific_enrichments.R


# Official analysis ends here
# ---------------------------



####################################################
# The above SNPs include all the GWAS in the paper;
# here we throw in a few of the latest brain GWAS
####################################################

# Note: You must run the above scripts first, and then 
# run these additional analyses.

# Munge summary statistics for bonus brain GWAS
bash scripts/gwas_sumstat.sh

# Perform LD score regression to partition heritability
# on a few new GWAS
bash scripts/cell_type_spe_regression_bonus_gwas.sh

# Make heritability enrichment plots like the one in the paper,
# for additional GWAS only
Rscript scripts/plot_cellspecific_enrichments_new_gwas.R

# Make the plots for additional GWAS and the ones included in the
# paper, combined
Rscript scripts/plot_cellspecific_enrichments_new_and_old_gwas.R
