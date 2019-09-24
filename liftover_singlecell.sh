# Transform single cell Atac-seq peaks bed files in from hg38 to hg19 

#Idr Optimal Peaks
atacseq_folder=$1
for x in astrocytes doublets excitatory_neurons inhibitory_neurons microglia neurons_unknown nigral_neurons oligodendrocytes opcs ; do
echo "Processing ${x}"
liftOver $atacseq_folder/$x.idr.optimal_peak.narrowPeak.gz  \
hg38ToHg19.over.chain.gz \
$atacseq_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed \
$atacseq_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 

bgzip $atacseq_folder/$x.idr.optimal_peak.narrowPeak.hg19.bed
bgzip $atacseq_folder/$x.idr.optimal_peak.narrowPeak.unmapp.bed 
done
