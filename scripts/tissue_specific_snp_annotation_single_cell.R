library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)

plink_dir = "data/1kg-plink-files"

read_bim = function(plink_dir){
	bim_fn_list=list.files(plink_dir,'bim',full.names=TRUE)
	bim=foreach(i=seq_along(bim_fn_list),.combine='rbind')%dopar%{
		bim_fn=bim_fn_list[i]
		message(bim_fn)
		fread(bim_fn,select=1:4,col.names=c('CHR','SNP','CM','BP'))
	}
	return(bim)
}

annotate = function(bim,peak,tissues){
	annot = bim
	setcolorder(annot,c('CHR','BP','SNP','CM'))
	annot[,c('start','end'):=list(BP,BP)]
	setkey(annot,CHR, start, end)

	for (t in tissues){
		peak_t = peak[str_detect(tissue,t),]
		setkey(peak_t,chr,start,end)
		overlap = foverlaps(annot,peak_t,by.x = c('CHR','start','end'), by.y = c('chr','start','end'), nomatch = 0)
		annot[,new := as.integer(SNP %in% overlap$SNP)]
		setnames(annot,'new',t)
	}
	annot$start = NULL
	annot$end = NULL
	return(annot)
}


############################################################
# Cell type specific annotations,
# first for IDR peaks and then for naive optimal peak calls
############################################################


tissues = c('astrocytes',  'excitatory_neurons', 'inhibitory_neurons', 'microglia', 'neurons_unknown', 'dopaminergic_neurons', 'oligodendrocytes', 'opcs')

#Idr Optimal Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificIdr/nodoublets/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.idr.optimal_peak.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}


# Overlap Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/CelltypeSpecificNaiveOverlapOptimalPeaksBedtoolsMerge/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificNaiveOverlap/nodoublets/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.overlap.optimal_peak.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}

############################################################
# Cell type clustered annotations,
# first for IDR peaks and then for naive optimal peak calls
############################################################

tissues = c('C10', 'C11', 'C12', 'NeuronC13', 'NeuronC14', 'NeuronC15', 'NeuronC17', 'NeuronC18', 'NeuronC19', 'NeuronC20', 'NeuronC21', 'NeuronC24', 'NeuronC25', 'NeuronC26', 'NeuronC27', 'NeuronC28', 'NeuronC29', 'NeuronC2', 'NeuronC30', 'NeuronC3', 'NeuronC6', 'NeuronC7', 'NeuronC8', 'NeuronC9')
#Idr Optimal Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/cluster_frags/idr_peaks/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/cluster_frags/idr_peaks/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.idr.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}


# Overlap Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/cluster_frags/overlap_peaks/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/cluster_frags/overlap_peaks/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.overlap.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}


############################################################
# Cell type subgroup annotations,
# first for IDR peaks and then for naive optimal peak calls
############################################################

tissues = c('LAMP5_Interneurons', 'Parvalbumin_Interneurons', 'Somatostatin_Interneurons', 'Striatal_Interneurons', 'Striatopallidal_MediumSpinyNeurons', 'vGLUT2neg_ExcitatoryNeurons', 'VIPIP_Interneurons', 'BDNFneg_ExcitatoryNeurons', 'BDNFpos_ExcitatoryNeurons', 'Striatonigral_MediumSpinyNeurons')

#Idr Optimal Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/group_frags/idr_peaks/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/group_frags/idr_peaks/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.idr.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}


# Overlap Peaks no doublets separate per tissue
peak_dir = 'data/atac-seq/group_frags/overlap_peaks/'
out_dir =  'output/ld_score_regression/tissue_specific_snp_annotation/group_frags/overlap_peaks/'
if (!dir.exists(out_dir)) {dir.create(out_dir)}
for (i in 1:length(tissues)){
	peak_fn = paste0(peak_dir, tissues[i],'.overlap.narrowPeak.hg19.bed.gz')
	peak = fread(peak_fn,select = 1:3, col.names = c('chr','start','end'))
	peak$tissue = tissues[i]
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues[i])

	for (j in 1:22){
		out = annot[CHR == j]
		out_fn = sprintf('%s/%s.chr%s.annot',out_dir,tissues[i],j)
		fwrite(out,out_fn,sep='\t')
	}
}


