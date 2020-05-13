##################################################
# Run all steps in LDSC regression for cell-type 
# specific heritability partitioning
##################################################

source scripts/load_local_env.sh

# Munge summary statistics
bash scripts/gwas_sumstat.sh

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

