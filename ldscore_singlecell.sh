export plink_dir=/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/
export annot_dir=/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/tissue_specific_snp_annotation/
export hapmap_dir=/srv/persistent/bliu2/shared/ldscore/hapmap3_snps/


export out_dir=/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/ldscore_merged/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/
export OPENBLAS_NUM_THREADS=15

for x in astrocytes excitatory_neurons inhibitory_neurons microglia neurons_unknown nigral_neurons oligodendrocytes opcs ;
do 	
		seq 1 22 | parallel -j 1  "python /users/bliu2/tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $annot_dir/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/$x.chr{}.annot --out $out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
done
