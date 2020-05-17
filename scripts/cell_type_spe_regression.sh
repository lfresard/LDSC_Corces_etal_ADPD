unset DISPLAY XAUTHORITY

export sumstats_folder="output/ld_score_regression/gwas_sumstats/"

# Load pre-munged GWAS sumstats files, in LDSC format
export parkinsons_sumstats_fn="$sumstats_folder/parkinsons_23andMe.sumstats.gz"
export alzheimers_sumstats_Kunkle_fn="$sumstats_folder/alzheimers_Kunkle.sumstats.gz"
export anorexia_fn="$sumstats_folder/AnorexiaNervosa_Duncan_2017.sumstats.gz"
export adhd_fn="$sumstats_folder/Attention_Deficit_2017.sumstats.gz"
export anxiety_fn="$sumstats_folder/Anxiety_Otowa_2016.sumstats.gz"
export neuroticism_fn="$sumstats_folder/Neuroticism_Symptoms_Okbay_2016.sumstats.gz"
export schizophrenia_fn="$sumstats_folder/schizophrenia_Li2017.sumstats.gz"
export leanbodymass_fn="$sumstats_folder/Lean_Body_Mass_Zillikens_2017.sumstats.gz"
export bone_density_fn="$sumstats_folder/Bone_Mineral_Density_Kemp_2017.sumstats.gz"
export cad_fn="$sumstats_folder/Coronary_Artery_Disease_Howson_2017.sumstats.gz"
export epilepsy_fn="$sumstats_folder/Epilepsy_Anney_2014.sumstats.gz"

export out_dir='output/ld_score_regression/partition_heritability_merged/'
export tissue_specific_annotation_dir='output/ld_score_regression/ldscore_merged'
export weight_dir="data/ld-score-weights"
export frq_dir="data/1kg-freqs"
export baseline_annotation_dir='data/ldsc-baseline'


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

python tools/ldsc/ldsc.py \
    --h2-cts $4 \
    --ref-ld-chr $baseline_annotation_dir/baseline. \
    --out $2/$3.merged \
    --ref-ld-chr-cts $5 \
    --w-ld-chr $weight_dir/weights.
}

export -f cell_type_SPE

for gwas in $parkinsons_sumstats_fn $alzheimers_sumstats_Kunkle_fn $anorexia_fn $adhd_fn $anxiety_fn $neuroticism_fn $schizophrenia_fn $leanbodymass_fn $bone_density_fn $cad_fn $epilepsy_fn
do
	# IDR peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificIdr/nodoublets/ \
	$out_dir/CelltypeSpecificIdr/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$out_dir/CelltypeSpecificIdr/nodoublets/CelltypeSpecificIdr.ldcts

	# Overlap peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	$out_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$out_dir/CelltypeSpecificNaiveOverlap/nodoublets/CelltypeSpecificOverlap.ldcts

	# IDR group peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/group_frags/idr_peaks/ \
	$out_dir/group_frags/idr_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/group_frags/idr_peaks/CelltypeSpecific.group.idr.ldcts

	# Overlap group peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/group_frags/overlap_peaks/ \
	$out_dir/group_frags/overlap_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/group_frags/overlap_peaks/CelltypeSpecific.group.overlap.ldcts

	# IDR cluster peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/cluster_frags/idr_peaks// \
	$out_dir/cluster_frags/idr_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/cluster_frags/idr_peaks/CelltypeSpecific.cluster.idr.ldcts

	# Overlap cluster peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/cluster_frags/overlap_peaks// \
	$out_dir/cluster_frags/overlap_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/cluster_frags/overlap_peaks/CelltypeSpecific.cluster.overlap.ldcts


done


new_gwas_dir='output/ld_score_regression/gwas_sumstats_new/'

alcohol_sumstats_fn=$new_gwas_dir/Alcohol-Dependence_Sanchez-Roige_2018.sumstats.gz
autism_sumstats_fn=$new_gwas_dir/Autism_Psychiatric-Genomics-Consortium_2017.sumstats.gz
bipolar_sumstats_fn=$new_gwas_dir/Bipolar-Disorder_Stahl_2019.sumstats.gz
depression_sumstats_fn=$new_gwas_dir/Depression_Howard_2019.sumstats.gz
epilepsy_new_sumstats_fn=$new_gwas_dir/Epilepsy_ILAE_2018.sumstats.gz
insomnia_sumstats_fn=$new_gwas_dir/Insomnia_Jansen_2019.sumstats.gz
ptsd_sumstats_fn=$new_gwas_dir/Post-Traumatic-Stress-Disorder_Duncan_2017.sumstats.gz
parkinsons_age_onset_sumstats_fn=$new_gwas_dir/Parkinsons-Age-At-Onset_Blauwendraat_2019.sumstats.gz
ocd_sumstats_fn=$new_gwas_dir/Obsessive-Compulsive-Disorder_Arnold_2017.sumstats.gz


for gwas in $alcohol_sumstats_fn $autism_sumstats_fn $bipolar_sumstats_fn $depression_sumstats_fn $epilepsy_new_sumstats_fn $insomnia_sumstats_fn $ptsd_sumstats_fn $parkinsons_age_onset_sumstats_fn $ocd_sumstats_fn
do
	# IDR peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificIdr/nodoublets/ \
	$out_dir/CelltypeSpecificIdr/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$out_dir/CelltypeSpecificIdr/nodoublets/CelltypeSpecificIdr.ldcts

	# Overlap peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	$out_dir/CelltypeSpecificNaiveOverlap/nodoublets/ \
	`basename $gwas` \
	$gwas \
	$out_dir/CelltypeSpecificNaiveOverlap/nodoublets/CelltypeSpecificOverlap.ldcts

	# IDR group peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/group_frags/idr_peaks/ \
	$out_dir/group_frags/idr_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/group_frags/idr_peaks/CelltypeSpecific.group.idr.ldcts

	# Overlap group peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/group_frags/overlap_peaks/ \
	$out_dir/group_frags/overlap_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/group_frags/overlap_peaks/CelltypeSpecific.group.overlap.ldcts

	# IDR cluster peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/cluster_frags/idr_peaks// \
	$out_dir/cluster_frags/idr_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/cluster_frags/idr_peaks/CelltypeSpecific.cluster.idr.ldcts

	# Overlap cluster peaks
	cell_type_SPE \
	$tissue_specific_annotation_dir/cluster_frags/overlap_peaks// \
	$out_dir/cluster_frags/overlap_peaks/ \
	`basename $gwas` \
	$gwas \
	$out_dir/cluster_frags/overlap_peaks/CelltypeSpecific.cluster.overlap.ldcts

done
