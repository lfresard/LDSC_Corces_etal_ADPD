# LD score regression: Cell type specific analyses

The full set of analyses for LD Score regression are outlined in the Bash script `run_all.sh`.

We have used relative paths wherever possible, but you will likely need to swap some of these
out to reference whatever paths your data are stored in (symlinks usually work well for this purpose).

GWAS data, pre-formatted for the analysis, are available at 
https://zenodo.org/record/3817811#.XrW9InVKhhE

---------------------

There are two ways of calling peaks : IDR and naive overlap

There are three peak groupings: basic cell type specific, cell type subgroups, and cell type clusters.

This makes a total of 2 x 3 = 6 analyses, so you will see each step in the analysis performed
six times. The results are not too different between the IDR and overlap peaks, so if you want quick
results, we recommend just running the steps for IDR peaks with basic cell-type-specific clustering.

## Process


A brief outline of what is done at each step:

### Format GWAS summary statistics for running LDSC regression

`bash scripts/gwas_sumstat.sh`

Summary statistics from the relevant GWAS are munged into a format
suitable for LDSC regression. If you download the GWAS files we linked on the
Zotero website, this step won't be necessary because they are already in
LDSC format (for build hg19).

### Extract the necessary columns from the ATAC-seq data

`bash scripts/prep_atac_seq.sh`

Subset ATAC-seq peaks BED files to the columns required for 
LD score analysis (only chromosome, start, end columns are needed,
for a total of three columns)

### Lift over ATAC-seq peaks from hg38 to hg19

`bash scripts/liftover_singlecell.sh`

Lift over ATAC-seq peaks from hg38 to hg19; obviously not necessary
if the peaks are already in hg19 form.

### Annotate SNPs by which tissues they belong to

`Rscript scripts/tissue_specific_snp_annotation_single_cell.R`

SNPs to be used in LD score regression are marked with all the tissues in which
they overlap with an ATAC-seq peak.

### Generate LD scores for each single cell type

`bash scripts/ldscore_singlecell.sh`

Generates LD scores for every single cell type, one at a time.

### Perform LD score regression to partition heritability

`bash scripts/cell_type_spe_regression.sh`

The actual step of partitioning heritability.

### Plot final enrichments

`Rscript scripts/plot_cellspecific_enrichments.R`

Making the plots shown in the paper.

## Required files and tools

### Tools

- Local installation of ldsc.py (should be placed in a `bin/` folder, and must be run using Python 2.7)
- UCSC liftOver tool, if ATAC-seq peaks not already in hg19 format

### Required LDSC regression files
- 1000 Genomes allele frequencies and PLINK-generated LD reference (see LDSC documentation at https://github.com/bulik/ldsc)
- Baseline LD score model for comparison, along with pre-computed LD-score weights (see LDSC documentation)
- Hapmap reference indicating which SNPs will be kept for LD Score regression (see LDSC documentation)
- UCSC chain files, if needing to use our liftOver script to convert ATAC-seq peaks from hg38 to hg19

### Data files for analysis
- ATAC-seq peaks in bed format
- GWAS in format produced by the LDSC `munge_sumstats.py` script, in the same
  genome build as the ATAC-seq peaks

