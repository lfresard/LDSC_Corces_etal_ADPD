res_files=list.files(pattern=".merged.cell_type_results.txt")

data_list=lapply(res_files, function(x){read.csv(file=x,header=T, sep="\t")})

gwas_list=c("adhd", "alzheimers", "anorexia", "anxiety", "bone_density", "cad", "epilepsy", "lean_body_mass","neuroticism", "parkinsons","skizophrenia")
data_list_gwas=Map(cbind, data_list, gwas = gwas_list)
library(tidyverse)

data.df=do.call(rbind,data_list_gwas)

data.df$gwas=factor(data.df$gwas, levels=c("lean_body_mass","cad","bone_density","adhd",  "neuroticism", "anorexia","epilepsy", "skizophrenia","anxiety","parkinsons","alzheimers"))
bigplot=ggplot(data.df, aes(x=gwas,y=-log10(Coefficient_P_value), color=Name)) +
	labs(y="-log10(Coefficient P-value)")+
	geom_point(size=3) + coord_flip()+ scale_color_brewer(palette="Set3")+ theme_bw()

ggsave("resuts_sc_gwas.pdf", bigplot,dpi = 300)
