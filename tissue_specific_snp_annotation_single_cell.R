library(data.table)
library(foreach)
library(doMC)
registerDoMC(20)
library(stringr)

source("local_paths.R")

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

if (FALSE)
{
	#Idr Optimal Peaks no doublets
	peak_fn = paste0(dir,'data/kundaje_version/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/allregions.idr.optimal_peak.narrowPeak.nodoublets.hg19.bed')
	tissues = c('astrocytes',  'excitatory_neurons', 'inhibitory_neurons', 'microglia', 'neurons_unknown', 'dopaminergic_neurons', 'oligodendrocytes', 'opcs')
	out_dir =  paste0(dir,'/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/nodoublets/')
	if (!dir.exists(out_dir)) {dir.create(out_dir)}

	peak = fread(peak_fn,select = 1:4, col.names = c('chr','start','end','tissue'))
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues)

	for (i in 1:22){
		out = annot[CHR == i]
		out_fn = sprintf('%s/chr%s.annot',out_dir,i)
		fwrite(out,out_fn,sep='\t')
	}
}


#Idr Optimal Peaks no doublets separate per tissue
peak_dir = paste0(dir,'data/kundaje_version/CelltypeSpecificIdrOptimalPeaksBedtoolsMerge/')
tissues = c('astrocytes',  'excitatory_neurons', 'inhibitory_neurons', 'microglia', 'neurons_unknown', 'dopaminergic_neurons', 'oligodendrocytes', 'opcs')
out_dir =  paste0(dir,'/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificIdr/nodoublets/')
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

if (FALSE)
{
	# Overlap Peaks no doublets
	peak_fn = paste0(dir,'data/kundaje_version/CelltypeSpecificNaiveOverlapOptimalPeaksBedtoolsMerge/allregions.overlap.optimal_peak.narrowPeak.hg19.bed')
	tissues = c('astrocytes',  'excitatory_neurons', 'inhibitory_neurons', 'microglia', 'neurons_unknown', 'dopaminergic_neurons', 'oligodendrocytes', 'opcs')
	out_dir =  paste0(dir,'/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificNaiveOverlap/nodoublets/')
	if (!dir.exists(out_dir)) {dir.create(out_dir)}

	peak = fread(peak_fn,select = 1:4, col.names = c('chr','start','end','tissue'))
	peak[,chr:=as.integer(str_replace(chr,'chr',''))]
	bim = read_bim(plink_dir)

	annot = annotate(bim,peak,tissues)

	for (i in 1:22){
		out = annot[CHR == i]
		out_fn = sprintf('%s/chr%s.annot',out_dir,i)
		fwrite(out,out_fn,sep='\t')
	}
}


# Overlap Peaks no doublets separate per tissue
peak_dir = paste0(dir,'data/kundaje_version/CelltypeSpecificNaiveOverlapOptimalPeaksBedtoolsMerge/')
tissues = c('astrocytes',  'excitatory_neurons', 'inhibitory_neurons', 'microglia', 'neurons_unknown', 'dopaminergic_neurons', 'oligodendrocytes', 'opcs')
out_dir =  paste0(dir,'/processed_data/ld_score_regression/tissue_specific_snp_annotation/CelltypeSpecificNaiveOverlap/nodoublets/')
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


