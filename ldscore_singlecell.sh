ad_pd='/oak/stanford/groups/smontgom/mgloud/projects/ad-pd-ldsc/'

export plink_dir=/srv/persistent/bliu2/shared/ldscore/1000G_plinkfiles/
export overlap_annot_dir=$ad_pd/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificNaiveOverlap/nodoublets/
export idr_annot_dir=$ad_pd/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificIdr/nodoublets
export hapmap_dir=/srv/persistent/bliu2/shared/ldscore/hapmap3_snps/


export overlap_out_dir=$ad_pd/processed_data/ld_score_regression/ldscore_merged/CelltypeSpecificNaiveOverlap/nodoublets/
export idr_out_dir=$ad_pd/processed_data/ld_score_regression/ldscore_merged/CelltypeSpecificIdr/nodoublets/
export OPENBLAS_NUM_THREADS=15

for x in astrocytes excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ;
do 	
	seq 1 22 | parallel -j 1  "python /users/bliu2/tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $overlap_annot_dir/$x.chr{}.annot --out $overlap_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
	seq 1 22 | parallel -j 1  "python /users/bliu2/tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $idr_annot_dir/$x.chr{}.annot --out $idr_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
done
