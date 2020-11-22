# Transforms GWAS summary statistics into right format for LDSC
# This one includes a few bonus GWAS traits that are quite new
# and were not included in the main analysis

new_gwas_dir='data/gwas/new-gwas'

hapmap_fn='data/munging/hapmap/w_hm3.snplist'

out_dir='output/ld_score_regression/gwas_sumstats_new/'

mkdir -p $out_dir

source activate python_ldsc

# Try a few new ones, just for kicks
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Alcohol-Dependence_Sanchez-Roige_2018/Alcohol-AUDIT-Total-Score.txt.gz --out $out_dir/Alcohol-Dependence_Sanchez-Roige_2018 --merge-alleles $hapmap_fn --N 121604 --a1 effect_allele --a2 non_effect_allele
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Autism_Psychiatric-Genomics-Consortium_2017/Autism-Worldwide.txt.gz  --out $out_dir/Autism_Psychiatric-Genomics-Consortium_2017 --merge-alleles $hapmap_fn --N-cas 16539 --N-con 157234 --a1 effect_allele --a2 non_effect_allele --ignore or
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Bipolar-Disorder_Stahl_2019/Bipolar-Disorder_Stahl_2019.txt.gz --out $out_dir/Bipolar-Disorder_Stahl_2019 --merge-alleles $hapmap_fn --N-cas 185285 --N-con 439741  --a1 effect_allele --a2 non_effect_allele
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Depression_Howard_2019/Depression-Meta-Analysis.txt.gz --out $out_dir/Depression_Howard_2019 --merge-alleles $hapmap_fn --N-cas 246363 --N-con 561190 --a1 effect_allele --a2 non_effect_allele --ignore or
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Epilepsy_ILAE_2018/All-Epilepsy.txt.gz --out $out_dir/Epilepsy_ILAE_2018 --merge-alleles $hapmap_fn --N-cas 15212 --N-con 29677 --a1 effect_allele --a2 non_effect_allele --ignore or
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Insomnia_Jansen_2019/Insomnia.txt.gz --out $out_dir/Insomnia_Jansen_2019 --merge-alleles $hapmap_fn --N 1331010 --a1 effect_allele --a2 non_effect_allele
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Post-Traumatic-Stress-Disorder_Duncan_2017/Post-Traumatic-Stress-Disorder_Duncan_2017.txt.gz --out $out_dir/Post-Traumatic-Stress-Disorder_Duncan_2017 --merge-alleles $hapmap_fn --N-cas 5182 --N-con 15548 --a1 effect_allele --a2 non_effect_allele --ignore or
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Parkinsons-Age-At-Onset_Blauwendraat_2019/Parkinsons-Age-At-Onset_Blauwendraat_2019.txt.gz --out $out_dir/Parkinsons-Age-At-Onset_Blauwendraat_2019 --merge-alleles $hapmap_fn --N 28568 --a1 effect_allele --a2 non_effect_allele --ignore MarkerName
python2.7 tools/ldsc/munge_sumstats.py --sumstats $new_gwas_dir/Obsessive-Compulsive-Disorder_Arnold_2017/Obsessive-Compulsive-Disorder_Arnold_2017.txt.gz --out $out_dir/Obsessive-Compulsive-Disorder_Arnold_2017 --merge-alleles $hapmap_fn --N-cas 2688 --N-con 7037 --a1 effect_allele --a2 non_effect_allele


