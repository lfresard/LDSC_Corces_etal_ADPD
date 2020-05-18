
# Load data shared across all peak groups
export plink_dir='data/1kg-plink-files'
export hapmap_dir='data/hapmap'

export OPENBLAS_NUM_THREADS=15

############################################################
# Cell type specific annotations,
# first for naive optimal peak calls and then for IDR
############################################################

# Load intermediate pre-processed files
export overlap_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificNaiveOverlap/nodoublets/'
export idr_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificIdr/nodoublets'

# Output directory locations
export overlap_out_dir='output/ld_score_regression/ldscore_merged/CelltypeSpecificNaiveOverlap/nodoublets/'
export idr_out_dir='output/ld_score_regression/ldscore_merged/CelltypeSpecificIdr/nodoublets/'

for x in astrocytes excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ;
do 	
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $overlap_annot_dir/$x.chr{}.annot --out $overlap_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $idr_annot_dir/$x.chr{}.annot --out $idr_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
done

############################################################
# Cell type subgroups annotations,
# first for naive optimal peak calls and then for IDR
############################################################

# Load intermediate pre-processed files
export group_frags_overlap_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/group_frags/overlap_peaks/'
export group_frags_idr_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/group_frags/idr_peaks/'

# Output directory locations
export group_frags_overlap_out_dir='output/ld_score_regression/ldscore_merged/group_frags/overlap_peaks/'
export group_frags_idr_out_dir='output/ld_score_regression/ldscore_merged/group_frags/idr_peaks/'

for x in LAMP5_Interneurons Parvalbumin_Interneurons Somatostatin_Interneurons Striatal_Interneurons Striatopallidal_MediumSpinyNeurons vGLUT2neg_ExcitatoryNeurons VIPIP_Interneurons BDNFneg_ExcitatoryNeurons BDNFpos_ExcitatoryNeurons Striatonigral_MediumSpinyNeurons;
do 	
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $group_frags_overlap_annot_dir/$x.chr{}.annot --out $group_frags_overlap_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $group_frags_idr_annot_dir/$x.chr{}.annot --out $group_frags_idr_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
done

############################################################
# Cell type clustered annotations,
# first for naive optimal peak calls and then for IDR
############################################################

# Load intermediate pre-processed files
export cluster_frags_overlap_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/cluster_frags/overlap_peaks/'
export cluster_frags_idr_annot_dir='output/ld_score_regression/tissue_specific_snp_annotation/cluster_frags/idr_peaks/'

# Output directory locations
export cluster_frags_overlap_out_dir='output/ld_score_regression/ldscore_merged/cluster_frags/overlap_peaks/'
export cluster_frags_idr_out_dir='output/ld_score_regression/ldscore_merged/cluster_frags/idr_peaks/'

for x in C10 C11 C12 NeuronC13 NeuronC14 NeuronC15 NeuronC17 NeuronC18 NeuronC19 NeuronC20 NeuronC21 NeuronC24 NeuronC25 NeuronC26 NeuronC27 NeuronC28 NeuronC29 NeuronC2 NeuronC30 NeuronC3 NeuronC6 NeuronC7 NeuronC8 NeuronC9;
for x in C11;
do 	
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $cluster_frags_overlap_annot_dir/$x.chr{}.annot --out $cluster_frags_overlap_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
	seq 1 22 | parallel -j 1  "python tools/ldsc/ldsc.py --l2 --bfile $plink_dir/1000G.mac5eur.{} --ld-wind-cm 1 --annot $cluster_frags_idr_annot_dir/$x.chr{}.annot --out $cluster_frags_idr_out_dir/$x.chr{} --print-snps $hapmap_dir/hm.{}.snp"
done
