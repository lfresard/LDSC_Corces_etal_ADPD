unset DISPLAY XAUTHORITY

#export alzheimers_sumstats_Lambert_fn='/srv/persistent/bliu2/neuro_finemap/processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/alzheimers_Lambert.sumstats.gz'
#export alzheimers_sumstats_Jansen_fn='/srv/persistent/bliu2/neuro_finemap/processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/alzheimers_Jansen.sumstats.gz'
export parkinsons_sumstats_fn='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/parkinsons_23andMe.sumstats.gz'
export alzheimers_sumstats_Kunkle_fn='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/alzheimers_Kunkle.sumstats.gz'
export anorexia_fn='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/AnorexiaNervosa_Duncan_2017.sumstats.gz'
export adhd_fn='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/Attention_Deficit_2017.sumstats.gz'
export anxiety_fn='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats//Anxiety_Otowa_2016.sumstats.gz'
export neuroticism_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats//Neuroticism_Symptoms_Okbay_2016.sumstats.gz"
export schizophrenia_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/schizophrenia_Li2017.sumstats.gz"
export leanbodymass_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/Lean_Body_Mass_Zillikens_2017.sumstats.gz"
export bone_density_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/Bone_Mineral_Density_Kemp_2017.sumstats.gz"
export cad_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/Coronary_Artery_Disease_Howson_2017.sumstats.gz"
export epilepsy_fn="/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats//Epilepsy_Anney_2014.sumstats.gz"

export out_dir='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/'
export tissue_specific_annotation_dir='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/ldscore_merged'
export weight_dir='/srv/persistent/bliu2/shared/ldscore/weights_hm3_no_hla/'
export frq_dir='/srv/persistent/bliu2/shared/ldscore/1000G_frq/'
export baseline_annotation_dir='/srv/persistent/bliu2/shared/ldscore/baseline/'


cell_type_SPE(){
mkdir -p $2
echo INFO - in_dir: $1
echo INFO - out_dir: $2
echo INFO - trait: $3
echo INFO - gwas: $4
echo INFO - ldcts: $5
if [[ -e "$2/$3.merged.results" ]]; then

echo "INFO - Already finished. Skipping."
return 0

fi 

python /users/bliu2/tools/ldsc/ldsc.py \
    --h2-cts $4 \
    --ref-ld-chr $baseline_annotation_dir/baseline. \
    --out $2/$3.merged \
    --ref-ld-chr-cts $5 \
    --w-ld-chr $weight_dir/weights.
}

export -f cell_type_SPE


cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
alzheimers_Kunkle_sc \
$alzheimers_sumstats_Kunkle_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
anorexia_Duncan_2017_sc \
$anorexia_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
parkinsons_sc \
$parkinsons_sumstats_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
adhd_sc \
$adhd_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts


cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
anxiety_sc \
$anxiety_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
neuroticism_sc \
$neuroticism_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
schizophrenia_Li2017 \
$schizophrenia_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
Lean_Body_Mass_Zillikens_2017 \
$leanbodymass_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
Bone_Mineral_Density_Kemp_2017 \
$bone_density_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
Coronary_Artery_Disease_Howson_2017 \
$cad_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

cell_type_SPE \
$tissue_specific_annotation_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
$out_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/ \
Epilepsy_Anney_2014 \
$epilepsy_fn \
/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/CelltypeSpecificIdr.ldcts

