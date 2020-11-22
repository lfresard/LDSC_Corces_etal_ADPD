library(tidyverse)

# This is ugly visually, but it makes sure the matching is done correctly
filenames = c()
traitnames = c()

filenames = c(filenames, "Alcohol-Dependence_Sanchez-Roige_2018.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "alcoholism")

filenames = c(filenames, "Autism_Psychiatric-Genomics-Consortium_2017.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "autism")

filenames = c(filenames, "Bipolar-Disorder_Stahl_2019.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "bipolar")

filenames = c(filenames, "Depression_Howard_2019.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "depression")

filenames = c(filenames, "Insomnia_Jansen_2019.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "insomnia")

filenames = c(filenames, "Obsessive-Compulsive-Disorder_Arnold_2017.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "ocd")

filenames = c(filenames, "Post-Traumatic-Stress-Disorder_Duncan_2017.sumstats.gz.merged.cell_type_results.txt")
traitnames = c(traitnames, "ptsd")

# Folders containining partitioned heritability results
ph_folders = c("output/ld_score_regression/partition_heritability_merged/group_frags/idr_peaks/",
	       "output/ld_score_regression/partition_heritability_merged/CelltypeSpecificIdr/nodoublets/",
	       "output/ld_score_regression/partition_heritability_merged/cluster_frags/idr_peaks/",
	       "output/ld_score_regression/partition_heritability_merged/group_frags/overlap_peaks/",
	       "output/ld_score_regression/partition_heritability_merged/cluster_frags/overlap_peaks/",
	       "output/ld_score_regression/partition_heritability_merged/CelltypeSpecificNaiveOverlap/nodoublets/")

set_names = c("group_idr",
	      "original_idr",
	      "cluster_idr",
	      "group_overlap",
	      "cluster_overlap",
	      "original_overlap")

for (i in 1:length(ph_folders))
{
	group_folder = ph_folders[i]
	group = set_names[i]

	full_filenames = paste0(group_folder, filenames)

	data_list=lapply(full_filenames, function(x){read.csv(file=x,header=T, sep="\t")})

	data_list_gwas=Map(cbind, data_list, gwas = traitnames)

	data.df=do.call(rbind,data_list_gwas)

	data.df$point_size = ifelse(data.df$Coefficient_P_value < 0.05, 3, 1)

	data.df$gwas=factor(data.df$gwas, levels=c("alcoholism", "autism", "bipolar", "depression", "insomnia", "ocd", "ptsd"))
	bigplot=ggplot(data.df, aes(x=gwas,y=-log10(Coefficient_P_value), color=Name, size=point_size)) +
		labs(y="-log10(Coefficient P-value)")+
		geom_hline(yintercept=-log10(0.05), linetype="dotted", color="red", size=1) +
		scale_size(range=c(2,5), guide = 'none') +
		geom_point() + coord_flip()+ 
		#scale_color_brewer(palette="Set3")+ 
		theme_bw()

	ggsave(paste0("output/ld_score_regression/plots/results_sc_new_gwas_", group, ".pdf"), bigplot,dpi = 300)

}


