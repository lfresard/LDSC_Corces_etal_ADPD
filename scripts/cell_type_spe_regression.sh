unset DISPLAY XAUTHORITY

export ad_pd='/oak/stanford/groups/smontgom/mgloud/projects/ad-pd-ldsc'

#export alzheimers_sumstats_Lambert_fn='/srv/persistent/bliu2/neuro_finemap/processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/alzheimers_Lambert.sumstats.gz'
#export alzheimers_sumstats_Jansen_fn='/srv/persistent/bliu2/neuro_finemap/processed_data/disease_enrichment/ld_score_regression/gwas_sumstats/alzheimers_Jansen.sumstats.gz'
export parkinsons_sumstats_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/parkinsons_23andMe.sumstats.gz"
export alzheimers_sumstats_Kunkle_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/alzheimers_Kunkle.sumstats.gz"
export anorexia_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/AnorexiaNervosa_Duncan_2017.sumstats.gz"
export adhd_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Attention_Deficit_2017.sumstats.gz"
export anxiety_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Anxiety_Otowa_2016.sumstats.gz"
export neuroticism_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Neuroticism_Symptoms_Okbay_2016.sumstats.gz"
export schizophrenia_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/schizophrenia_Li2017.sumstats.gz"
export leanbodymass_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Lean_Body_Mass_Zillikens_2017.sumstats.gz"
export bone_density_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Bone_Mineral_Density_Kemp_2017.sumstats.gz"
export cad_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Coronary_Artery_Disease_Howson_2017.sumstats.gz"
export epilepsy_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Epilepsy_Anney_2014.sumstats.gz"

export out_dir=$ad_pd/processed_data/ld_score_regression/partition_heritability_merged/
export tissue_specific_annotation_dir=$ad_pd/processed_data/ld_score_regression/ldscore_merged
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

export parkinsons_sumstats_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/parkinsons_23andMe.sumstats.gz"
export alzheimers_sumstats_Kunkle_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/alzheimers_Kunkle.sumstats.gz"
export anorexia_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/AnorexiaNervosa_Duncan_2017.sumstats.gz"
export adhd_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Attention_Deficit_2017.sumstats.gz"
export anxiety_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Anxiety_Otowa_2016.sumstats.gz"
export neuroticism_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Neuroticism_Symptoms_Okbay_2016.sumstats.gz"
export schizophrenia_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/schizophrenia_Li2017.sumstats.gz"
export leanbodymass_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Lean_Body_Mass_Zillikens_2017.sumstats.gz"
export bone_density_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Bone_Mineral_Density_Kemp_2017.sumstats.gz"
export cad_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats/Coronary_Artery_Disease_Howson_2017.sumstats.gz"
export epilepsy_fn="$ad_pd/processed_data/ld_score_regression/gwas_sumstats//Epilepsy_Anney_2014.sumstats.gz"



for gwas in $parkinsons_sumstats_fn $alzheimers_sumstats_Kunkle_fn $anorexia_fn $adhd_fn $anxiety_fn $neuroticism_fn $schizophrenia_fn $leanbodymass_fn $bone_density_fn $cad_fn $epilepsy_fn
do
	# IDR peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificIdr/nodoublets/ \
	$out_dir/CelltypeSpecificIdr/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$ad_pd/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdr/nodoublets/CelltypeSpecificIdr.ldcts

	# Overlap peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	$out_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$ad_pd/processed_data/ld_score_regression/partition_heritability_merged/CelltypeSpecificNaiveOverlap/nodoublets/CelltypeSpecificOverlap.ldcts


done

