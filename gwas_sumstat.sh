# Transforms GWAS summary statistics into right format for LDSC

out_dir='/users/lfresard/scratch/alzheimers_parkinsons/processed_data/ld_score_regression/gwas_sumstats/'
hapmap_fn='/srv/persistent/bliu2/shared/ldscore/w_hm3.snplist'

alzheimers_gwas_Kunkle_fn='/users/lfresard/scratch/alzheimers_parkinsons/data/gwas/alzheimers/Kunkle2019/Kunkle_etal_Stage1_results.txt.gz'
parkinsons_23andMe='/users/lfresard/scratch/alzheimers_parkinsons/data/gwas/parkinsons/23andMe/GWAS_sumstat.txt'
schizophrenia_Li2017='/mnt/lab_data/montgomery/shared/gwas/Schizophrenia_Li_2017/Schizophrenia-Li-2017.txt.gz' # problem in file
AnorexiaNervosa_Duncan_2017='/mnt/lab_data/montgomery/shared/gwas/Anorexia-Nervosa_Duncan_2017/pgc.ed.freeze1.summarystatistics.July2017.txt'
Attention_Deficit_2017='/mnt/lab_data/montgomery/shared/gwas/Attention-Deficit-Hyperactivity-Disorder_Demontis_2017/adhd_jul2017.gz'
Anxiety_Otowa_2016='/mnt/lab_data/montgomery/shared/gwas/Anxiety_Otowa_2016/anxiety.meta.full.fs.tbl'
Neuroticism_Symptoms_Okbay_2016='/mnt/lab_data/montgomery/shared/gwas/Neuroticism-Symptoms_Okbay_2016/Neuroticism_Full.txt.gz'
Lean_Body_Mass_Zillikens_2017="/mnt/lab_data/montgomery/shared/gwas/Lean-Body-Mass_Zillikens_2017/wholebodyleanmass.results.metal_.txt.gz"
Bone_Mineral_Density_Kemp_2017="/mnt/lab_data/montgomery/shared/gwas/Bone-Mineral-Density-Estimated_Kemp_2017/BEurope-Bmd-As-C-Gwas-SumStats.txt_0.gz"
Coronary_Artery_Disease_Howson_2017="/mnt/lab_data/montgomery/shared/gwas/Coronary-Artery-Disease_Howson_2017/Howson-JMM_CHD_Mixed_2017.txt.gz"
Epilepsy_Anney_2014="/mnt/lab_data/montgomery/shared/gwas/Epilepsy_Anney_2014/ILAE_All_Epi_11.8.14.txt.gz"

mkdir -p $out_dir



source activate python_ldsc

/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $alzheimers_gwas_Kunkle_fn --out $out_dir/alzheimers_Kunkle --merge-alleles $hapmap_fn --N-cas 21982 --N-con 41944
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $parkinsons_23andMe --out $out_dir/parkinsons_23andMe --merge-alleles $hapmap_fn --N-cas 6476 --N-con 302042
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $schizophrenia_Li2017 --out $out_dir/schizophrenia_Li2017 --merge-alleles $hapmap_fn --N-cas 34241  --N-con 45604
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $AnorexiaNervosa_Duncan_2017 --out $out_dir/AnorexiaNervosa_Duncan_2017 --merge-alleles $hapmap_fn --N-cas 3495  --N-con 10982
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Attention_Deficit_2017 --out $out_dir/Attention_Deficit_2017 --merge-alleles $hapmap_fn --N-cas 20183  --N-con 35191
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Anxiety_Otowa_2016 --out $out_dir/Anxiety_Otowa_2016 --merge-alleles $hapmap_fn --N-col totalN
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Neuroticism_Symptoms_Okbay_2016 --out $out_dir/Neuroticism_Symptoms_Okbay_2016 --merge-alleles $hapmap_fn --N 170911
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Lean_Body_Mass_Zillikens_2017 --out $out_dir/Lean_Body_Mass_Zillikens_2017 --merge-alleles $hapmap_fn --N 38292
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Bone_Mineral_Density_Kemp_2017 --out $out_dir/Bone_Mineral_Density_Kemp_2017 --merge-alleles $hapmap_fn --N 142487 --a1 ALLELE0 --a2 ALLELE1
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Coronary_Artery_Disease_Howson_2017 --out $out_dir/Coronary_Artery_Disease_Howson_2017 --merge-alleles $hapmap_fn --N-cas 88192 --N-con 162544 --a1 effect_allele --a2 other_allele
/users/bliu2/tools/ldsc/munge_sumstats.py --sumstats $Epilepsy_Anney_2014 --out $out_dir/Epilepsy_Anney_2014 --merge-alleles $hapmap_fn --N-cas 8696 --N-con 26157 --a1 Allele1 --a2 Allele2


