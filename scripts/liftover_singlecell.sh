# Transform single cell Atac-seq peaks bed files in from hg38 to hg19

# For this process to work, you must have liftOver installed and in your $PATH;
# you must also have the appropriate chain file present at data/chain-files/hg38ToHg19.over.chain.gz

############################################################
# Cell type specific analysis,
# first for IDR peaks and then for naive optimal peak calls
############################################################

export atacseq_idr_folder="data/atac-seq/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/"
export atacseq_overlap_folder="data/atac-seq/CelltypeSpecificNaiveOverlapOptimalPeaksBedtoolsMerge/"

# Idr Optimal Peaks
for x in astrocytes doublets excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ; do
	echo "Processing ${x}"
	liftOver $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed \
	$atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 

	bgzip $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed
	bgzip $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 
done


# Naive overlap peaks
for x in astrocytes doublets excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ; do
	echo "Processing ${x}"
	liftOver $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.hg19.bed \
	$atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.unmapp.bed 

	bgzip $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.hg19.bed
	bgzip $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.unmapp.bed 
done

############################################################
# Cell type subgroup analysis,
# first for IDR peaks and then for naive optimal peak calls
############################################################

export atacseq_group_frags_idr_folder="data/atac-seq/group_frags/idr_peaks"
export atacseq_group_frags_overlap_folder="data/atac-seq/group_frags/overlap_peaks"

#Idr Optimal Peaks
for x in LAMP5_Interneurons Parvalbumin_Interneurons Somatostatin_Interneurons Striatal_Interneurons Striatopallidal_MediumSpinyNeurons vGLUT2neg_ExcitatoryNeurons VIPIP_Interneurons BDNFneg_ExcitatoryNeurons BDNFpos_ExcitatoryNeurons Striatonigral_MediumSpinyNeurons; do
	echo "Processing ${x}"
	liftOver $atacseq_group_frags_idr_folder/$x.idr.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_group_frags_idr_folder/$x.idr.narrowPeak.hg19.bed \
	$atacseq_group_frags_idr_folder/$x.idr.narrowPeak.unmapp.bed 

	bgzip $atacseq_group_frags_idr_folder/$x.idr.narrowPeak.hg19.bed
	bgzip $atacseq_group_frags_idr_folder/$x.idr.narrowPeak.unmapp.bed 
done

# Naive overlap peaks
for x in LAMP5_Interneurons Parvalbumin_Interneurons Somatostatin_Interneurons Striatal_Interneurons Striatopallidal_MediumSpinyNeurons vGLUT2neg_ExcitatoryNeurons VIPIP_Interneurons BDNFneg_ExcitatoryNeurons BDNFpos_ExcitatoryNeurons Striatonigral_MediumSpinyNeurons; do
	echo "Processing ${x}"
	liftOver $atacseq_group_frags_overlap_folder/$x.overlap.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_group_frags_overlap_folder/$x.overlap.narrowPeak.hg19.bed \
	$atacseq_group_frags_overlap_folder/$x.overlap.narrowPeak.unmapp.bed 

	bgzip $atacseq_group_frags_overlap_folder/$x.overlap.narrowPeak.hg19.bed
	bgzip $atacseq_group_frags_overlap_folder/$x.overlap.narrowPeak.unmapp.bed 
done


############################################################
# Cell type cluster analysis,
# first for IDR peaks and then for naive optimal peak calls
############################################################


export atacseq_cluster_frags_idr_folder="data/atac-seq/cluster_frags/idr_peaks"
export atacseq_cluster_frags_overlap_folder="data/atac-seq/cluster_frags/overlap_peaks"


#Idr Optimal Peaks
for x in C10 C11 C12 NeuronC13 NeuronC14 NeuronC15 NeuronC17 NeuronC18 NeuronC19 NeuronC20 NeuronC21 NeuronC24 NeuronC25 NeuronC26 NeuronC27 NeuronC28 NeuronC29 NeuronC2 NeuronC30 NeuronC3 NeuronC6 NeuronC7 NeuronC8 NeuronC9; do
	echo "Processing ${x}"
	liftOver $atacseq_cluster_frags_idr_folder/Cluster-$x.idr.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_cluster_frags_idr_folder/$x.idr.narrowPeak.hg19.bed \
	$atacseq_cluster_frags_idr_folder/$x.idr.narrowPeak.unmapp.bed 

	bgzip -f $atacseq_cluster_frags_idr_folder/$x.idr.narrowPeak.hg19.bed
	bgzip -f $atacseq_cluster_frags_idr_folder/$x.idr.narrowPeak.unmapp.bed 
done


# Naive overlap peaks
for x in C10 C11 C12 NeuronC13 NeuronC14 NeuronC15 NeuronC17 NeuronC18 NeuronC19 NeuronC20 NeuronC21 NeuronC24 NeuronC25 NeuronC26 NeuronC27 NeuronC28 NeuronC29 NeuronC2 NeuronC30 NeuronC3 NeuronC6 NeuronC7 NeuronC8 NeuronC9; do
	echo "Processing ${x}"
	liftOver $atacseq_cluster_frags_overlap_folder/Cluster-$x.overlap.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_cluster_frags_overlap_folder/$x.overlap.narrowPeak.hg19.bed \
	$atacseq_cluster_frags_overlap_folder/$x.overlap.narrowPeak.unmapp.bed 

	bgzip -f $atacseq_cluster_frags_overlap_folder/$x.overlap.narrowPeak.hg19.bed
	bgzip -f $atacseq_cluster_frags_overlap_folder/$x.overlap.narrowPeak.unmapp.bed 
done


############################################################
# All clusters regardless of whether neuronal,
# first for IDR peaks and then for naive optimal peak calls
############################################################


export atacseq_cluster_complete_frags_idr_folder="data/atac-seq/ClusterSpecificIDROptimalPeaksBedtoolsMerge"
export atacseq_cluster_complete_frags_overlap_folder="data/atac-seq/ClusterSpecificNaiveOverlapOptimalPeaksBedtoolsMerge"


#Idr Optimal Peaks
for x in `seq 1 24`; do
	echo "Processing ${x}"
	liftOver $atacseq_cluster_complete_frags_idr_folder/Cluster$x.idr.col-subs.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_cluster_complete_frags_idr_folder/Cluster$x.idr.narrowPeak.hg19.bed \
	$atacseq_cluster_complete_frags_idr_folder/Cluster$x.idr.narrowPeak.unmapp.bed 

	bgzip -f $atacseq_cluster_complete_frags_idr_folder/Cluster$x.idr.narrowPeak.hg19.bed
	bgzip -f $atacseq_cluster_complete_frags_idr_folder/Cluster$x.idr.narrowPeak.unmapp.bed 
done


# Naive overlap peaks
for x in `seq 1 24`; do
	echo "Processing ${x}"
	liftOver $atacseq_cluster_complete_frags_overlap_folder/Cluster$x.optimal_peak.col-subs.narrowPeak.gz  \
	data/chain-files/hg38ToHg19.over.chain.gz	\
	$atacseq_cluster_complete_frags_overlap_folder/Cluster$x.overlap.narrowPeak.hg19.bed \
	$atacseq_cluster_complete_frags_overlap_folder/Cluster$x.overlap.narrowPeak.unmapp.bed 

	bgzip -f $atacseq_cluster_complete_frags_overlap_folder/Cluster$x.overlap.narrowPeak.hg19.bed
	bgzip -f $atacseq_cluster_complete_frags_overlap_folder/Cluster$x.overlap.narrowPeak.unmapp.bed 
done


