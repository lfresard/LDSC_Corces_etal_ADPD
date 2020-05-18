####################################################################
# This script pulls out the columns from the ATAC-seq files that
# will be required to annotate SNPs; all other columns are
# irrelevant to us for this analysis
#
# Note: this step may or may not actually be necessary depending
# on what format the data are obtained in. The important thing
# is that after this step you should have a BED-formatted file
# for each cell type
#
# where column 1 shows chromosome, column 2 shows peak start,
# and column 3 shows peak end
####################################################################



# Cluster level

for cluster in C10 C11 C12 NeuronC13 NeuronC14 NeuronC15 NeuronC17 NeuronC18 NeuronC19 NeuronC20 NeuronC21 NeuronC24 NeuronC25 NeuronC26 NeuronC27 NeuronC28 NeuronC29 NeuronC2 NeuronC30 NeuronC3 NeuronC6 NeuronC7 NeuronC8 NeuronC9; do

	echo $cluster

	# Format the IDR peaks files
	zcat data/atac-seq/raw/clusterFrags/idr_peaks/Cluster-$cluster.idr.narrowPeak.gz | cut -f1,2,3 | sort -k2,2n > data/atac-seq/cluster_frags/idr_peaks/Cluster-$cluster.idr.narrowPeak
	bgzip -f data/atac-seq/cluster_frags/idr_peaks/Cluster-$cluster.idr.narrowPeak

	# Format the overlap peaks files
	zcat data/atac-seq/raw/clusterFrags/overlap_peaks/Cluster-$cluster.overlap.narrowPeak.gz | cut -f1,2,3 | sort -k2,2n > data/atac-seq/cluster_frags/overlap_peaks/Cluster-$cluster.overlap.narrowPeak
	bgzip -f data/atac-seq/cluster_frags/overlap_peaks/Cluster-$cluster.overlap.narrowPeak
done


# Group level

for cluster in LAMP5_Interneurons Parvalbumin_Interneurons Somatostatin_Interneurons Striatal_Interneurons Striatopallidal_MediumSpinyNeurons vGLUT2neg_ExcitatoryNeurons VIPIP_Interneurons BDNFneg_ExcitatoryNeurons BDNFpos_ExcitatoryNeurons Striatonigral_MediumSpinyNeurons; do

	echo $cluster

	# Format the IDR peaks files
	zcat data/atac-seq/raw/groupFrags/idr_peaks/$cluster.idr.narrowPeak.gz | cut -f1,2,3 | sort -k2,2n > data/atac-seq/group_frags/idr_peaks/$cluster.idr.narrowPeak
	bgzip -f data/atac-seq/group_frags/idr_peaks/$cluster.idr.narrowPeak

	# Format the overlap peaks files
	zcat data/atac-seq/raw/groupFrags/overlap_peaks/$cluster.overlap.narrowPeak.gz | cut -f1,2,3 | sort -k2,2n > data/atac-seq/group_frags/overlap_peaks/$cluster.overlap.narrowPeak
	bgzip -f data/atac-seq/group_frags/overlap_peaks/$cluster.overlap.narrowPeak
done
