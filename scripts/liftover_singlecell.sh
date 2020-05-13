# Transform single cell Atac-seq peaks bed files in from hg38 to hg19 

export atacseq_idr_folder=$ad_pd/data/kundaje_version/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/
export atacseq_overlap_folder=$ad_pd/data/kundaje_version/CelltypeSpecificNaiveOverlapOptimalPeaksBedtoolsMerge/

#Idr Optimal Peaks
for x in astrocytes doublets excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ; do
	echo "Processing ${x}"
	liftOver $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.gz  \
	$chain_file	\
	$atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed \
	$atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 

	bgzip $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed
	bgzip $atacseq_idr_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 
done


# Naive overlap peaks
for x in astrocytes doublets excitatory_neurons inhibitory_neurons microglia neurons_unknown dopaminergic_neurons oligodendrocytes opcs ; do
	echo "Processing ${x}"
	liftOver $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.gz  \
	$chain_file	\
	$atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.hg19.bed \
	$atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.unmapp.bed 

	bgzip $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.hg19.bed
	bgzip $atacseq_overlap_folder/$x.overlap.optimal_peak.narrowPeak.unmapp.bed 
done
